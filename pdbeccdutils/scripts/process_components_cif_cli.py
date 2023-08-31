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
"""
Script for PDBeChem backend infrastructure.
Processes the wwPDB Chemical Components Dictionary file `components.cif`
producing files for:

http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/

To do this components.cif is split into individual PDB chemical component
definitions cif files, sdf files, pdb files and image files.
In addition creates chem_comp.xml and chem_comp.list for all components.
More detailed description can be found here:

https://gitlab.ebi.ac.uk/pdbe/release/pdbechem
"""
import argparse
import logging
import os
from pathlib import Path

import pdbeccdutils
import rdkit
from pdbeccdutils.core import ccd_reader, ccd_writer, prd_reader, prd_writer
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.depictions import DepictionManager
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.helpers.cif_tools import (
    get_ccd_cif_enriched_dir,
    get_prd_cif_enriched_dir,
    get_prdcc_code,
)
from pdbeccdutils.core.fragment_library import FragmentLibrary
from pdbeccdutils.core.models import ConformerType, DepictionSource
from pdbeccdutils.utils import config
from pdbeccdutils.utils.pubchem_downloader import PubChemDownloader


def is_valid_path(path):
    p = Path(path)

    if not p.exists():
        raise argparse.ArgumentTypeError(f"{path} does not exists.")

    return p


class PDBeChemManager:
    """Manager orchestrating computation and generation of all parts of the
    PDBeChem update process.
    """

    def __init__(
        self,
        pubchem_templates="",
        general_templates=config.general_templates,
        library_path=config.fragment_library,
        procedure="ccd",
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

        self.procedure = procedure

    def run(self, components_path, out_dir):
        """Process components

        Args:
            components_path (Path): Path to the components.cif file.
            out_dir (Path): Path to the out_dir
        """
        logging.info("Reading in component")
        if self.procedure == "ccd":
            reader_result = ccd_reader.read_pdb_cif_file(components_path)
            component = reader_result.component
            component_output_dir = get_ccd_cif_enriched_dir(component.id, out_dir)
            outfile_prefix = component.id

        else:
            reader_result = prd_reader.read_pdb_cif_file(components_path)
            component = reader_result.component
            component_output_dir = get_prd_cif_enriched_dir(component.id, out_dir)
            outfile_prefix = get_prdcc_code(component.id)

        os.makedirs(component_output_dir, exist_ok=True)
        self.process_single_component(
            reader_result, component_output_dir, outfile_prefix
        )
        logging.info(f"Processing of {component.id} is complete")

    def process_single_component(self, ccd_reader_result, out_dir, outfile_prefix):
        """Process single PDB-CCD component.

        Args:
            ccd_reader_result (CCDReaderResult): pdbeccdutils parser output.
            out_dir (Path): Out directory
            outfile_prefix: Prefix of output filename

        Return:
            bool: Whether or not all the files were succesfully written.
        """
        component = ccd_reader_result.component
        logging.info(f"{component.id} | processing...")

        # check parsing and conformer degeneration
        self._check_component_parsing(ccd_reader_result)
        self._generate_ideal_structure(component)

        # download templates if the user wants them.
        if self.pubchem is not None:
            self._download_template(component)

        # search fragment library
        self._search_fragment_library(component)

        # get scaffolds
        self._compute_component_scaffolds(component)

        # write out files
        self._generate_depictions(component, out_dir, outfile_prefix)
        self._export_structure_formats(component, out_dir, outfile_prefix)

    def _check_component_parsing(self, ccd_reader_result):
        """Checks components parsing and highlights issues encountered with
        the molecule: errors/warnings during the parsing process,
        unrecoverable sanitization issues, inchikey mismatch between what
        was in the source file and is reproduced by rdkit.

        Args:
            ccd_reader_result (CCDReaderResult): Output of the parsing process.
        """

        if ccd_reader_result.warnings:
            logging.debug(f'warnings: {";".join(ccd_reader_result.warnings)}')

        if ccd_reader_result.errors:
            logging.debug(f'errors: {";".join(ccd_reader_result.errors)}')

        if not ccd_reader_result.sanitized:
            logging.debug("sanitization issue.")

        if not ccd_reader_result.component.inchikey_from_rdkit_matches_ccd():
            logging.debug("inchikey mismatch.")

    def _download_template(self, component: Component):
        """Attempts to download a pubchem template for the given component

        Args:
            component (Component): Component to be used.
        """
        component_downloaded = self.pubchem.process_template(component)
        if component_downloaded:
            logging.debug("downloaded new pubchem template.")

    def _generate_ideal_structure(self, component: Component):
        """Checks whether or not the component has degenerated ideal
        coordinates. If so, new conformer is attempted to be generated.

        Args:
            component (Component): Component to be
                processed.
        Return:
            bool: Whether the ideal coordinates have been successfully
            recalculated, false otherwise.
        """
        result = component.compute_3d()

        if component.has_degenerated_conformer(ConformerType.Ideal):
            logging.debug(f"{component.id} has degenerated ideal coordinates.")

        if not result:
            logging.debug(f"{component.id} has error in generating 3D conformation.")

        return result

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

    def _generate_depictions(
        self, component: Component, out_dir: str, outfile_prefix: str
    ):
        """Generate nice 2D depictions for the component and
        depiction annotations in JSON format. Presently depictions
        are generated in the following resolutions (100,200,300,400,500)
        with and without atom names.

        Args:
            component (Component): Component to be depicted.
            out_dir (str): Where the depictions should be stored.
            outfile_prefix (str): Prefix of output filename
        """
        depiction_result = component.compute_2d(self.depictions)

        if depiction_result.source == DepictionSource.Failed:
            logging.info("failed to generate 2D image.")
        else:
            if depiction_result.score > 0.99:
                logging.info("collision free image could not be generated.")
            logging.info(
                f"2D generated using {depiction_result.source.name} with score {depiction_result.score}."
            )

        wedge_bonds = depiction_result.template_name != "cube"

        for i in range(100, 600, 100):
            component.export_2d_svg(
                os.path.join(out_dir, f"{outfile_prefix}_{i}.svg"),
                width=i,
                wedge_bonds=wedge_bonds,
            )
            component.export_2d_svg(
                os.path.join(out_dir, f"{outfile_prefix}_{i}_names.svg"),
                width=i,
                names=True,
                wedge_bonds=wedge_bonds,
            )

        component.export_2d_annotation(
            os.path.join(out_dir, f"{outfile_prefix}_annotation.json"),
            wedge_bonds=wedge_bonds,
        )

    def _export_structure_formats(
        self, component: Component, out_dir: Path, outfile_prefix: str
    ):
        """Writes out component in a different formats as required for the
        PDBeChem FTP area.

        Args:
            component (Component): Component being processed.
            out_dir (Path): Where the results should be written
            outfile_prefix (str): Prefix of output filename
        """

        self.__write_molecule(
            out_dir / f"{outfile_prefix}_model.sdf",
            component,
            False,
            ConformerType.Model,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}_ideal.sdf",
            component,
            False,
            ConformerType.Ideal,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}_ideal_alt.pdb",
            component,
            True,
            ConformerType.Ideal,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}_model_alt.pdb",
            component,
            True,
            ConformerType.Model,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}_ideal.pdb",
            component,
            False,
            ConformerType.Ideal,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}_model.pdb",
            component,
            False,
            ConformerType.Model,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}.cml",
            component,
            False,
            ConformerType.Model,
        )
        self.__write_molecule(
            out_dir / f"{outfile_prefix}.cif",
            component,
            False,
            ConformerType.AllConformers,
        )

    def __write_molecule(self, path, component, alt_names, conformer_type):
        """Write out deemed structure.

        Args:
            path (Path): Path where the molecule will be stored.
            component (Component): Component to be written.
            alt_names (bool): Whether or not molecule will be written with
                alternate names.
            conformer_type (Component): Conformer to be written.
        """
        try:
            if self.procedure == "ccd":
                ccd_writer.write_molecule(
                    path,
                    component,
                    remove_hs=False,
                    alt_names=alt_names,
                    conf_type=conformer_type,
                )
            else:
                prd_writer.write_molecule(
                    path,
                    component,
                    remove_hs=False,
                    alt_names=alt_names,
                    conf_type=conformer_type,
                )
        except Exception:
            logging.error(f"error writing {path}.")

            with open(path, "w") as f:
                f.write("")


# region pre-light tasks
def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         argparse.Namespace parser
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-i", "--input-cif", help="Input PDB-CCD/PDB-PRD cif file")
    parser.add_argument(
        "-g",
        "--general-templates",
        type=is_valid_path,
        default=config.general_templates,
        help="Use general templates in SDF format instead of those supplied with the code.",
    )
    parser.add_argument(
        "-p",
        "--pubchem-templates",
        type=is_valid_path,
        help="Path to the directory with pubchem templates in sdf format.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=is_valid_path,
        required=True,
        help="Path to output directory",
    )
    parser.add_argument(
        "-fl",
        "--fragment-library",
        type=is_valid_path,
        default=config.fragment_library,
        help="Use this fragment library in place of the one supplied with the code.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Turn on debug message logging output"
    )

    # region parse procedure
    procedure = parser.add_mutually_exclusive_group()
    procedure.add_argument(
        "--ccd",
        action="store_const",
        const="ccd",
        dest="procedure",
        help="Pipeline is going to process CCD mmCIF files.",
    )
    procedure.add_argument(
        "--prd",
        action="store_const",
        const="prd",
        dest="procedure",
        help="Pipeline is going to process PRD mmCIF files.",
    )

    procedure.set_defaults(procedure="ccd")
    # endregion parse procedure

    return parser


# endregion


def main():
    """Runs the PDBeChem pipeline"""
    parser = create_parser()
    args = parser.parse_args()

    frm = "[%(asctime)-15s]  %(message)s"
    lvl = logging.DEBUG if __name__ == "__main__" and args.debug else logging.INFO

    logging.basicConfig(level=lvl, format=frm, datefmt="%a, %d %b %Y %H:%M:%S")
    logging.info("PDBeChem pipeline using:")
    logging.info(
        f"pdbeccdutils core v. {pdbeccdutils.__version__}, RDKit v. {rdkit.__version__}"
    )

    logging.info("Settings:")
    for k, v in vars(args).items():
        logging.info(f'{"":5s}{k:25s}{v}')

    pdbechem = PDBeChemManager(
        args.pubchem_templates,
        args.general_templates,
        args.fragment_library,
        args.procedure,
    )
    pdbechem.run(args.input_cif, args.output_dir)


if __name__ == "__main__":
    main()
