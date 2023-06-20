#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
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

"""Generic helper functions that may be re-used
"""

import os
import sys
import argparse
import logging
import requests
from requests.adapters import HTTPAdapter
from urllib3 import Retry

import rdkit

import pdbeccdutils


def bond_pdb_order(value_order):
    """
    Transpils mmcif bond order into rdkit language

    Args:
        value_order (str): bond type as a str

    Returns:
        rdkit.Chem.rdchem.BondType: -- bond type
    """
    if value_order.casefold() == "sing":
        return rdkit.Chem.rdchem.BondType(1)
    if value_order.casefold() == "doub":
        return rdkit.Chem.rdchem.BondType(2)
    if value_order.casefold() == "trip":
        return rdkit.Chem.rdchem.BondType(3)

    return None


def find_atom_index(mol, residue_id, atom_id):
    for atom in mol.GetAtoms():
        if (
            atom.GetProp("residue_id") == residue_id
            and atom.GetProp("component_atom_id") == atom_id
        ):
            return atom.GetIdx()


def get_additional_fields(auth_asym_id: str) -> tuple[str, str]:
    """Gets original auth_asym_id and assembly operator from auth_asym_id
    in assembly file

    Args:
        auth_asym_id: auth_asym_id value form mmCIF file

    Returns:
        A tuple of (orig_auth_asym_id, operator), where orig_auth_asym_id is the auth_asym_id
        in input structure and operator is the operator used by model-server

    """
    orig_auth_asym_id = auth_asym_id
    operator = ""
    if "_" in auth_asym_id:
        orig_auth_asym_id = auth_asym_id.split("_")[0]
        operator = "_" + auth_asym_id.split("_")[1]

    return (orig_auth_asym_id, operator)


def set_up_logger(args: argparse.Namespace) -> logging.Logger:
    """Sets up application level logging

    Args:
        args: Parsed arguments.

    Returns:
        Application logger.
    """
    level = logging.DEBUG if args.debug else logging.INFO
    frm = "%(asctime)-15s: %(module)s: %(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=frm, datefmt="%Y-%m-%d  %H:%M:%S")

    if not args.no_header:
        logging.info("PDBeBm pipeline using:")
        logging.info(
            f"pdbeccdutils core v. {pdbeccdutils.__version__}, RDKit v. {rdkit.__version__}"
        )
        logging.info("Settings:")
        for k, v in vars(args).items():
            logging.info(f"  {k:25s}{v}")


# parameter validation
def validate_path_exists(parameter: str, value: str) -> None:
    """Checks if the parameter exists.

    Args:
        parameter: Name of the input parameter
        value: Value of the input parameter
    """
    if not os.path.exists(value):
        print(
            f"Supplied parameter {parameter} is an invalid path {value}.",
            file=sys.stderr,
        )
        sys.exit(os.EX_NOINPUT)


def check_args(args: argparse.Namespace) -> None:
    """Validates supplied arguments.

    Args:
        args: An argparse namespace containing the required arguments
    """

    input_fields = ["input_cif", "general_templates"]
    args_values = vars(args)
    for field in input_fields:
        validate_path_exists(field, args_values[field])

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.input)

    if not args.pdb_id:
        args.pdb_id = os.path.basename(args.input)[0:4]


def requests_retry_session(
    retries=4,
    backoff_factor=1,
    status_forcelist=(429, 500, 502, 503, 504),
    session=None,
):

    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )

    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session
