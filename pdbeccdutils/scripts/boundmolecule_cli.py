# software from PDBe: Protein Data Bank in Europe; http://pdbe.org
#
# Copyright 2018 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#

import argparse
import logging
import os
import json
from pdbeccdutils.core import clc_reader, clc_writer
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.depictions import DepictionManager
from pdbeccdutils.core.exceptions import EntryFailedException
from pdbeccdutils.core.models import ConformerType, DepictionSource
from pdbeccdutils.core.fragment_library import FragmentLibrary
from pdbeccdutils.utils import config
from pdbeccdutils.helpers import cif_tools, helper
from pdbeccdutils.utils.pubchem_downloader import PubChemDownloader
from pdbeccdutils.core.boundmolecule import infer_bound_molecules


class PDBeBmManager:
    """Pipeline to identify bound-molecules and process"""

    def __init__(
        self,
        pubchem_templates=str(),
        general_templates=config.general_templates,
        library_path=config.fragment_library,
        discarded_ligands=config.DISCARDED_RESIDUES,
    ):
        """Initialize manager

        Args:
            pubchem_templates (str): Path to the pubchem templates
            general_templates (str, optional): Path to the general templates.
                Defaults to config.general_templates.
            library_path (str, optional): Path to the fragments library:
                Defaults to config.fragment_library.
        """
        # helper class to download templates if needed
        self.pubchem = (
            PubChemDownloader(pubchem_templates) if pubchem_templates else None
        )

        # helper class to get nice depictions
        self.depictions = DepictionManager(pubchem_templates, general_templates)

        # Fragments library to get substructure matches
        self.fragment_library = FragmentLibrary(library_path)

        self.discarded = discarded_ligands

    def process_entry(self, input_cif: str, pdb_id: str, output_dir: str):
        """Process a single entry in the pipeline:

        Args:
            pdb_id (str): PDB id
            output_dir: Path to output directory
        """

        bm_file = os.path.join(output_dir, "bound_molecules.json")
        logging.info(f"Processing {pdb_id}")
        logging.info(f"Working directory: {output_dir}")
        logging.info(f"Bound molecules file: {bm_file}")

        os.makedirs(output_dir, exist_ok=True)

        fixed_mmcif_file = os.path.join(output_dir, f"{pdb_id}_processed.cif")
        cif_tools.fix_updated_mmcif(input_cif, fixed_mmcif_file)
        self.process_boundmolecules(pdb_id, output_dir)

    def process_boundmolecules(self, pdb_id: str, output_dir: str):
        """Identifies and processes bound-molecules from a protein entry

        Args:
            pdb_id: PDB id
            output_dir: Path to output directory
        """
        fixed_mmcif_file = os.path.join(output_dir, f"{pdb_id}_processed.cif")
        if os.path.isfile(fixed_mmcif_file):
            clc_reader_results = []
            bms = infer_bound_molecules(
                fixed_mmcif_file, self.discarded, assembly=False
            )
            for i, bm in enumerate(bms, start=1):
                bm_id = f"bm{i}"
                reader_result = clc_reader.infer_multiple_chem_comp(
                    fixed_mmcif_file, bm, bm_id, sanitize=True
                )
                if reader_result:
                    clc_reader_results.append(reader_result)

            for i, clc_reader_result in enumerate(clc_reader_results, start=1):
                clc_id = f"CLC_{i}"
                component = clc_reader_result.component
                component.id = clc_id
                clc_out_dir = os.path.join(output_dir, component.id)
                os.makedirs(clc_out_dir, exist_ok=True)
                self.process_single_component(clc_reader_result, clc_out_dir)

            self._write_out_bm(pdb_id, bms, clc_reader_results, output_dir)
        else:
            raise EntryFailedException(f"Preprocessing of {pdb_id} failed")

    def process_single_component(
        self, clc_reader_result: clc_reader.CLCReaderResult, output_dir: str
    ):
        """Processes identified bound-molecules
            * Checks components parsing and highlights issues encountered with the molecule
            * Generates 3D conformer coordinates
            * Generated 2D depictions

        Args:
            clc_reader_result: List of CLCReaderResult
            output_dir: Path to ooutput directory
        """
        component = clc_reader_result.component
        logging.info(f"{component.id} | processing...")

        # check parsing
        self._check_component_parsing(clc_reader_result)
        self._generate_ideal_structure(component)

        # download templates if the user wants them.
        if self.pubchem is not None:
            self._download_template(component)

        # search fragment library
        self._search_fragment_library(component)

        # get scaffolds
        self._compute_component_scaffolds(component)

        # write out files
        self._generate_depictions(component, output_dir)
        self._export_structure_formats(component, output_dir)

    def _write_out_bm(
        self,
        pdb_id: str,
        bound_molecules,
        clc_reader_results: list[clc_reader.CLCReaderResult],
        output_dir: str,
    ):
        """Writes out the details of bound-molecules

        Args:
            pdb_id: PDB ID
            clc_reader_result: List of CLCReaderResult
            output_dir: Path to ooutput directory
        """
        result_bag = {
            "entry": pdb_id,
            "boundMolecules": [],
        }

        for i, bm in enumerate(bound_molecules, start=1):
            bm_id = f"bm{i}"
            inchi = ""
            inchikey = ""
            for clc_reader_result in clc_reader_results:
                component = clc_reader_result.component
                if bm.is_equivalent(clc_reader_result.bound_molecule):
                    inchi = component.inchi
                    inchikey = component.inchikey
            result_bag["boundMolecules"].append(
                {
                    "id": bm_id,
                    "composition": bm.to_dict(),
                    "inchi": inchi,
                    "inchikey": inchikey,
                }
            )
        bm_file = os.path.join(output_dir, "bound_molecules.json")
        with open(bm_file, "w") as f:
            json.dump(result_bag, f, sort_keys=True, indent=4)

    def _check_component_parsing(
        self, clc_reader_result: list[clc_reader.CLCReaderResult]
    ):
        """Checks components parsing and highlights issues encountered with
        the molecule: errors/warnings during the parsing process,
        unrecoverable sanitization issues

        Args:
            ccd_reader_result (CCDReaderResult): Output of the parsing process.
        """
        component = clc_reader_result.component
        for warning in clc_reader_result.warnings:
            logging.warning(f"{component.id} | {warning}")

        for error in clc_reader_result.errors:
            logging.error(f"{component.id} | {error}")

        if not clc_reader_result.sanitized:
            logging.warning(f"{component.id} | sanitization issue")

    def _download_template(self, component: Component):
        """Attempts to download a pubchem template for the given component

        Args:
            component (Component): Component to be used.
        """
        component_downloaded = self.pubchem.process_template(component)
        if component_downloaded:
            logging.debug(f"{component.id} | downloaded new pubchem template.")

    def _generate_ideal_structure(self, component: Component):
        """Generates 3D conformer coordinates and checks if the molecule has
        degenerated coordinates

        Args:
            component (Component): Component to be
                processed.
        Return:
            bool: Whether the ideal coordinates have been successfully
            recalculated, false otherwise.
        """

        result = component.compute_3d()
        if result:
            if component.has_degenerated_conformer(ConformerType.Computed):
                logging.debug(f"{component.id} has degenerated Computed coordinates.")

        if component.has_degenerated_conformer(ConformerType.Model):
            logging.debug(f"{component.id} has degenerated Model coordinates.")

        if not result:
            logging.debug(f"{component.id} has error in generating 3D conformation.")

        return result

    def _generate_depictions(self, component: Component, out_dir: str):
        """Generate nice 2D depictions for the component and
        depiction annotations in JSON format. Presently depictions
        are generated in the following resolutions (100,200,300,400,500)
        with and without atom names.

        Args:
            component (Component): Component to be depicted.
            out_dir (str): Where the depictions should be stored.
        """
        depiction_result = component.compute_2d(self.depictions)

        if depiction_result.source == DepictionSource.Failed:
            logging.debug(f"{component.id} failed to generate 2D image")
        else:
            if depiction_result.score > 0.99:
                logging.debug(
                    f"{component.id} collision free image could not be generated"
                )
            logging.debug(
                f"{component.id} 2D generated using {depiction_result.source.name} with score {depiction_result.score}."
            )

        wedge_bonds = depiction_result.template_name != "cube"

        for i in range(100, 600, 100):
            component.export_2d_svg(
                os.path.join(out_dir, f"{component.id}_{i}.svg"),
                width=i,
                wedge_bonds=wedge_bonds,
            )
            component.export_2d_svg(
                os.path.join(out_dir, f"{component.id}_{i}_names.svg"),
                width=i,
                names=True,
                wedge_bonds=wedge_bonds,
            )

        component.export_2d_annotation(
            os.path.join(out_dir, f"{component.id}_annotation.json"),
            wedge_bonds=wedge_bonds,
        )

    def _search_fragment_library(self, component: Component):
        """Search fragment library to find hits

        Args:
            component (Component): Component to be processed
        """

        matches = component.library_search(self.fragment_library)

        if matches:
            logging.debug(
                f"{len(matches)} matches found in the library `{self.fragment_library.name}`."
            )

    def _compute_component_scaffolds(self, component: Component):
        """Compute scaffolds for a given component.

        Args:
            component (Component): Component to be processed
        """

        try:
            component.get_scaffolds()
        except CCDUtilsError as e:
            logging.error(str(e))

            return

        logging.debug(f"{len(component.scaffolds)} scaffold(s) were found.")

    def _export_structure_formats(self, component: Component, out_dir: str):
        """Writes out component in a different formats as required for the
        PDBeChem FTP area.

        Args:
            component (Component): Component being processed.
            out_dir (Path): Where the results should be written
        """

        self.__write_molecule(
            os.path.join(out_dir, f"{component.id}_model.sdf"),
            component,
            ConformerType.Model,
        )

        self.__write_molecule(
            os.path.join(out_dir, f"{component.id}_model.pdb"),
            component,
            ConformerType.Model,
        )
        self.__write_molecule(
            os.path.join(out_dir, f"{component.id}.cml"),
            component,
            ConformerType.Model,
        )
        self.__write_molecule(
            os.path.join(out_dir, f"{component.id}.cif"),
            component,
            ConformerType.AllConformers,
        )

    def __write_molecule(self, path, component, conformer_type):
        """Write out deemed structure.

        Args:
            path (Path): Path where the molecule will be stored.
            component (Component): Component to be written.
            alt_names (bool): Whether or not molecule will be written with
                alternate names.
            conformer_type (Component): Conformer to be written.
        """
        try:
            clc_writer.write_molecule(
                path,
                component,
                remove_hs=False,
                conf_type=conformer_type,
            )
        except Exception:
            logging.error(f"error writing {path}.")

            with open(path, "w") as f:
                f.write("")


def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         argparse.Namespace parser
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-i", "--input-cif", required=True, help="Input structure in mmCIF format"
    )

    parser.add_argument(
        "-p",
        "--pubchem-templates",
        help="Path to the directory with pubchem templates in sdf format.",
    )
    parser.add_argument(
        "-g",
        "--general-templates",
        default=config.general_templates,
        help="Use general templates in SDF format instead of those supplied with the code.",
    )

    parser.add_argument(
        "-fl",
        "--fragment-library",
        default=config.fragment_library,
        help="Use this fragment library in place of the one supplied with the code.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        help="Create an output directory with files suitable for PDBeChem ftp directory",
    )
    parser.add_argument("-id", "--pdb-id", type=str, help="PDB ID", required=False)
    parser.add_argument(
        "--discarded-ligands",
        type=str,
        default=",".join(config.DISCARDED_RESIDUES),
        help="Comma separated name of ligands to be discarded",
        required=False,
    )

    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Turn off header information for the script.",
    )

    parser.add_argument(
        "--debug", action="store_true", help="Turn on debug message logging output"
    )

    return parser


def run():
    """Main method to execute pdbebm pipeline"""
    parser = create_parser()
    args = parser.parse_args()
    helper.check_args(args)
    discarded_ligands = config.DISCARDED_RESIDUES.union(
        {item.strip() for item in args.discarded_ligands.split(",")}
    )
    helper.set_up_logger(args)
    bpm = PDBeBmManager(
        args.pubchem_templates,
        args.general_templates,
        args.fragment_library,
        discarded_ligands,
    )
    bpm.process_entry(args.input_cif, args.pdb_id, args.output_dir)
