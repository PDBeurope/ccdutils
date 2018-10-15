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

from mmCif.mmcifIO import MMCIF2Dict

from arpeggio.core import InteractionComplex


class BoundMoleculeContainer:
    """Class to house bound molecules stored in the mmCIF file.
    Essentially it represents a graph, which can contain a multiple
    of disconnected components.
    """

    def __init__(self, residues=set(), connections=set()):
        self.residues = residues
        self.connections = connections

    def __str__(self):
        return '-'.join(map(lambda l: str(l), self.residues))

    def to_arpeggio_format(self):
        return list(map(lambda l: l.to_arpeggio(), self.residues))

    def add_connection(self, e):
        self.connections.add(e)

    def add_node(self, n):
        self.residues.add(n)

    def to_dict(self):
        nodes = sorted(list(self.residues), key=lambda l: (int(l.res_id), l.chain))

        results_bag = {'residues': [], 'connections': []}
        results_bag['residues'] = list(map(lambda l: l.to_dict(), nodes))
        results_bag['connections'] = list(map(lambda l: (nodes.index(l[0]), nodes.index(l[1])), self.connections))

        return results_bag

    def pop_bound_molecule(self):
        """Get bound molecule.

        Returns:
            BoundMoleculeContainer: connected component representing a
            single bound molecule.
        """
        if len(self.residues) == 0:
            return None

        visited_residues = set()
        visited_connections = set()
        stack = set([self.residues.pop()])

        while len(stack) > 0:
            processed_node = stack.pop()
            visited_residues.add(processed_node)
            edges = list(filter(lambda l: processed_node in l.residues, self.connections))

            for e in edges:
                visited_connections.add(e)
                self.connections.remove(e)
                other_node = e.get_other(processed_node)
                stack.add(other_node)
                self.residues.remove(other_node)

        return BoundMoleculeContainer(visited_residues, visited_connections)


class Residue:
    """Represents a single residue.
    """

    def __init__(self, name, chain, res_id):
        self.name = name
        self.chain = chain
        self.res_id = res_id

    def __eq__(self, other):
        if self.name != other.name:
            return False
        if self.chain != other.chain:
            return False
        if self.res_id != other.res_id:
            return False

        return True

    def to_dict(self):
        """Returns a dictionary representation of a given residue.

        Returns:
            (:obj:`dict` of :obj:`str`): Dictionary representation along
            with the mmCIF keys.
        """
        return {
            'label_comp_id': self.name,
            'auth_asym_id': self.chain,
            'auth_seq_id': self.res_id
        }

    def __hash__(self):
        return hash(self.chain + self.res_id + self.name)

    def __str__(self):
        return f'/{self.name}/{self.res_id}/{self.chain}/'

    def to_arpeggio(self):
        """Gets Arpeggio style representation of a residue e.g. `/A/129/`

        Returns:
            str: Residue description in Arpeggio style.
        """
        return f'/{self.chain}/{self.res_id}/'


class Connection:
    """Class representing a bond parsed fromt the _struct_conn namespace
    """

    def __init__(self, a, b):
        self.residues = (a, b)

    def __eq__(self, other):
        return self[0] in other.residues and self[1] in other.residues

    def __hash__(self):
        return hash(self[0]) * hash(self[1])

    def get_other(self, a):
        """Get the other residue from the bond.

        Args:
            a (Residue): A residue whose partner we are looking for

        Returns:
            Residue: The other residue which is a part of the connection.
        """
        if a in self.residues:
            return self[1] if self[0] == a else self[0]
        else:
            return None

    def __str__(self):
        return f'{self[0]} - {self[1]}'

    def __getitem__(self, i):
        """Indexer method to access n-th atom in the bond.

        Args:
            i (int): atom index

        Raises:
            AttributeError: if index is other than 0,1

        Returns:
            Residue: residue
        """
        if i < 0 or i > 1:
            raise AttributeError('Connection has just two partners.')
        return self.residues[i]


class ProtLigInteractions:
    def __init__(self, structure, to_discard=[]):
        """Create protein - ligand interaction object.

        Args:
            structure (str): Path to the structure to be processed
            to_discard (list, optional): Defaults to []. List of residue
                names to be discarded prior to protein-ligand interaction
                lookup.
        """
        self.bound_molecules = self._infer_bound_molecules(structure, to_discard)
        self.path = structure

    def get_all_interactions(self, interaction_cutoff=5.0, compensation_factor=0.1,
                             include_neighbours=False):
        """Retrieve interactions for all bound molecules found in the
        entry.

        Args:
           interaction_cutoff (float, optional): Defaults to 5.0. Distance
                cutoff for grid points to be `interacting` with the entity.
            compensation_factor (float, optional): Defaults to 0.1.
                Compensation factor for VdW radii dependent interaction types.
            include_neighbours (bool, optional): Defaults to False. Include
                non-bonding interactions between residues that are next to
                each other in sequence.

        Returns:
            list of dict of str: List of interactions in a dictionary like schema.
        """
        i = 0
        results = {}

        for bm in self.bound_molecules:
            i += 1
            selection = bm.to_arpeggio_format()

            results[f'bm{i}'] = self.get_interaction(selection, interaction_cutoff,
                                                     compensation_factor, include_neighbours)

        return results

    def get_interaction(self, selection, interaction_cutoff=5.0, compensation_factor=0.1,
                        include_neighbours=False):
        """Retrieve interactions in the protein within a given selection.

        Args:
            selection (str): Selection in the Arpeggio format. e.g. /A/123/
            interaction_cutoff (float, optional): Defaults to 5.0. Distance
                cutoff for grid points to be `interacting` with the entity.
            compensation_factor (float, optional): Defaults to 0.1.
                Compensation factor for VdW radii dependent interaction types.
            include_neighbours (bool, optional): Defaults to False. Include
                non-bonding interactions between residues that are next to
                each other in sequence.

        Returns:
            dict of str: Interactions in a dictionary like schema.
        """

        compl = InteractionComplex(self.path)
        compl.structure_checks()
        compl.address_ambiguities()
        compl.run_arpeggio(selection, interaction_cutoff, compensation_factor, include_neighbours)

        return compl.get_contacts()

    def _infer_bound_molecules(self, structure, to_discard):
        """Identify bound molecules in the input protein structure.

        Args:
            structure (str): Path to the structure.
            to_discard (list of str): List of residue names to be discarded

        Returns:
            [list of BoundMoleculeContainer]: Bound molecules found in the pdb entry.
        """
        bms = []
        temp = self.__parse_bound_molecules(structure, to_discard)

        while True:
            component = temp.pop_bound_molecule()
            if component is None:
                return bms
            else:
                bms.append(component)

        return bms

    def __parse_bound_molecules(self, path, to_discard):
        """Parse information from the information about HETATMS from the
        `_pdbx_nonpoly_scheme` (or `_atom_sites` if the former one is not
        present) and connectivity among them from `_struct_conn`.

        Args:
            path (str): Path to the mmCIF structure
            to_discard (list of str): List of residue names to be discarded.

        Returns:
            BoundMolecule: All the bound molecules in a given entry.
        """
        def __filter_ligands_from_nonpoly_schema(schema, bms):
            ligands = set([i for i in schema['mon_id']])
            bms.residues = set(filter(lambda l: l.name in ligands, bms.residues))

        def __parse_ligands_from_atom_sites(atom_sites):
            g = BoundMoleculeContainer()
            for i in range(len(atom_sites['id'])):
                if atom_sites['group_PDB'][i] == 'HETATM':
                    n = Residue(
                        atom_sites['label_comp_id'][i],  # aka label_comp_id
                        atom_sites['auth_asym_id'][i],  # aka auth_asym_id
                        atom_sites['auth_seq_id'][i])  # aka auth_seq_id

                    if n.name not in to_discard:
                        g.add_node(n)

            return g

        def __add_connections(g, struct_conn):
            for i in range(len(struct_conn['id'])):
                ptnr1 = filter(lambda l:
                               l.name == struct_conn['ptnr1_label_comp_id'][i] and
                               l.res_id == struct_conn['ptnr1_auth_seq_id'][i], g.residues)
                ptnr2 = filter(lambda l:
                               l.name == struct_conn['ptnr2_label_comp_id'][i] and
                               l.res_id == struct_conn['ptnr2_auth_seq_id'][i], g.residues)

                for x in ptnr1:
                    for y in ptnr2:
                        ptnr1_chain = struct_conn['ptnr1_auth_asym_id'][i]
                        ptnr2_chain = struct_conn['ptnr2_auth_asym_id'][i]

                        if x.chain == ptnr1_chain and y.chain == ptnr2_chain:
                            g.add_connection(Connection(x, y))
                            continue

                        if '-' in x.chain and '-' in y.chain:
                            x_split = x.chain.split('-')
                            y_split = y.chain.split('-')

                            if x_split[0] == ptnr1_chain and y_split[0] == ptnr2_chain and x_split[1] == y_split[1]:
                                g.add_connection(Connection(x, y))

        parsed_str = list(MMCIF2Dict().parse(path).values())[0]

        if '_pdbx_nonpoly_scheme' not in parsed_str:
            return BoundMoleculeContainer()

        bms = __parse_ligands_from_atom_sites(parsed_str['_atom_site'])
        __filter_ligands_from_nonpoly_schema(parsed_str['_pdbx_nonpoly_scheme'], bms)

        if '_struct_conn' in parsed_str:
            __add_connections(bms, parsed_str['_struct_conn'])

        return bms

# endregion
