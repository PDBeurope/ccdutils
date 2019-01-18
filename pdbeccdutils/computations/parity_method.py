"""Lightweight implementation of the parity method provided by
Jon Tyczak to Abhik to be used by the cofactors pipeline.
"""

from typing import NamedTuple

from rdkit import Chem
from rdkit.Chem import rdFMCS


class ParityResult(NamedTuple):
    """
    NamedTuple for the result of parity method along with the details
    necessary for calculating the similarity score.

    Attributes:
        template_atoms (int): Number of template atoms.
        query_atoms (int): Number of query molecule atoms.
        match_count (int): Size o the common subgaph match.
        similarity_score (float): Calculate similarity score.
"""
    template_atoms: int
    query_atoms: int
    match_count: int
    similarity_score: float


def _get_matches(mol, smarts):
    """Gets the subgraph match for the parity method.

    Args:
        mol(rdkit.Chem.rdchem.Mol): Molecule to be queried
        smarts(str): Molecule representation in SMARTS.

    Returns:
        int: Molecular subgraph matches.
    """

    patt = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(patt, uniquify=False)
    return matches


def _generate_sim_score(template, query, smarts):
    """Given the two moleculs and their common subgraph represented as
    smarts string returns similarity score.

    Args:
        template(rdkit.Chem.rdchem.Mol): template molecule.
        query(rdkit.Chem.Mol): query molecule.
        smarts(str): common subgraph in the SMARTS format.

    Returns:
        (int, int): tuple of size of the common subgraph and the
        similarity score.
    """

    if template.GetNumAtoms() == 1:
        smarts = Chem.MolToSmarts(template)
    elif query.GetNumAtoms() == 1:
        smarts = Chem.MolToSmarts(query)
    if smarts is None:
        best_matches = 0
        best_sim_score = 0.0
    else:
        matches_1 = _get_matches(template, smarts)
        matches_2 = _get_matches(query, smarts)
        best_matches = 0
        size_1 = template.GetNumAtoms()
        size_2 = query.GetNumAtoms()

        for match_1 in matches_1:
            for match_2 in matches_2:
                matches = 0
                for i in range(0, len(match_1)):
                    atom_1 = template.GetAtomWithIdx(match_1[i])
                    atom_2 = query.GetAtomWithIdx(match_2[i])
                    symbol_1 = atom_1.GetSymbol()
                    symbol_2 = atom_2.GetSymbol()
                    if symbol_1 == symbol_2:
                        matches += 1
                if matches > best_matches:
                    best_matches = matches
        best_sim_score = float(best_matches) / float(size_1 + size_2 - best_matches)

    return (best_matches, best_sim_score)


def compare_molecules(template, query, thresh=0.01):
    """Given the two molecules calculates their similarity score.
    If expected similarity score is lower than the threshold nothing
    is calculated.

    Args:
        template(rdkit.Chem.rdchem.Mol): Template molecule
        query(rdkit.Chem.rdchem.Mol): Query molecule
        thresh(float, optional): Defaults to 0.01: Threshold score for
            the match to be considered.

    Returns:
        ParityResult: Result of the PARITY comparison.
    """
    templ_atoms = template.GetNumAtoms()
    query_atoms = query.GetNumAtoms()

    min_num_atoms = min(templ_atoms, query_atoms)
    max_sim_score = float(min_num_atoms) / float(templ_atoms + query_atoms - min_num_atoms)

    if max_sim_score < thresh:
        return ParityResult(templ_atoms, query_atoms, 0, 0)

    mcs_graph = rdFMCS.FindMCS([template, query],
                               bondCompare=rdFMCS.BondCompare.CompareAny,
                               atomCompare=rdFMCS.AtomCompare.CompareAny,
                               timeout=40,
                               completeRingsOnly=True)

    matches, sim_score = _generate_sim_score(template, query, mcs_graph.smartsString)

    return ParityResult(templ_atoms, query_atoms, matches, sim_score)
