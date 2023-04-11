import argparse
import logging
import os
import json

from pdbeccdutils.core import bm_reader
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.depictions import DepictionManager
from pdbeccdutils.core.exceptions import EntryFailedException
from pdbeccdutils.core.models import ConformerType, DepictionSource
from pdbeccdutils.utils import config
from pdbeccdutils.helpers import cif_tools, helper
from pdbeccdutils.utils.pubchem_downloader import PubChemDownloader


class PDBeBmManager:
    """Manager orchestrating computation and generation of all parts of the
    PDBeBm update process.
    """

    def __init__(
        self,
        pubchem_templates=str(),
        general_templates=config.general_templates,
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
        self.discarded = discarded_ligands

    def process_entry(self, input_cif, pdb_id, output_dir):
        """Process a single entry in the pdbeboundmolecule pipeline:

        Args:
            pdb_id (str): PDB id
        """

        bm_file = os.path.join(output_dir, "bound_molecules.json")
        logging.info(f"Processing {pdb_id}")
        logging.info(f"Working directory: {output_dir}")
        logging.info(f"Bound molecules file: {bm_file}")

        os.makedirs(output_dir, exist_ok=True)

        fixed_mmcif_file = os.path.join(output_dir, f"{pdb_id}_processed.cif")
        cif_tools.fix_updated_mmcif(input_cif, fixed_mmcif_file)
        self.process_boundmolecules(pdb_id, output_dir)

    def process_boundmolecules(self, pdb_id, output_dir):
        fixed_mmcif_file = os.path.join(output_dir, f"{pdb_id}_processed.cif")
        if os.path.isfile(fixed_mmcif_file):
            bm_reader_results = bm_reader.read_pdb_updated_cif_file(
                fixed_mmcif_file, pdb_id, sanitize=True
            )
            for bm_reader_result in bm_reader_results:
                component = bm_reader_result.component
                bm_out_dir = os.path.join(output_dir, component.id)
                os.makedirs(bm_out_dir, exist_ok=True)
                self.process_bm_component(bm_reader_result, bm_out_dir)
            self._write_out_bm(pdb_id, bm_reader_results, output_dir)
        else:
            raise EntryFailedException(f"Preprocessing of {pdb_id} failed")

    def process_bm_component(self, bm_reader_result, output_dir):

        component = bm_reader_result.component
        logging.info(f"{component.id} | processing...")

        # check parsing
        self._check_component_parsing(bm_reader_result)
        self._generate_ideal_structure(component)

        # download templates if the user wants them.
        if self.pubchem is not None:
            self._download_template(component)

        self._generate_depictions(component, output_dir)

    def _write_out_bm(self, pdb_id, bm_reader_results, output_dir):
        result_bag = {
            "entry": pdb_id,
            "boundMolecules": [],
        }
        for bm_reader_result in bm_reader_results:
            component = bm_reader_result.component
            bm = bm_reader_result.bound_molecule
            result_bag["boundMolecules"].append(
                {
                    "id": component.id,
                    "composition": bm.to_dict(),
                    "inchi": component.inchi,
                    "inchikey": component.inchikey,
                }
            )

        bm_file = os.path.join(output_dir, "bound_molecules.json")
        with open(bm_file, "w") as f:
            json.dump(result_bag, f, sort_keys=True, indent=4)

    def _check_component_parsing(self, bm_reader_result):
        """Checks components parsing and highlights issues encountered with
        the molecule: errors/warnings during the parsing process,
        unrecoverable sanitization issues

        Args:
            ccd_reader_result (CCDReaderResult): Output of the parsing process.
        """
        component = bm_reader_result.component
        for warning in bm_reader_result.warnings:
            logging.warning(f"{component.id} | {warning}")

        for error in bm_reader_result.errors:
            logging.error(f"{component.id} | {error}")

        if not bm_reader_result.sanitized:
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
        if result:
            if component.has_degenerated_conformer(ConformerType.Computed):
                logging.debug("has degenerated Computed coordinates.")

        if component.has_degenerated_conformer(ConformerType.Model):
            logging.debug("has degenerated Model coordinates.")

        if not result:
            logging.debug("error in generating 3D conformation.")

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
            logging.debug(f"{component.id} | failed to generate 2D image")
        else:
            if depiction_result.score > 0.99:
                logging.debug(
                    f"{component.id} | collision free image could not be generated"
                )
            logging.debug(
                f"{component.id} | 2D generated using {depiction_result.source.name} with score {depiction_result.score}."
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
    bm_manager = PDBeBmManager(
        args.pubchem_templates, args.general_templates, discarded_ligands
    )
    bm_manager.process_entry(args.input_cif, args.pdb_id, args.output_dir)
