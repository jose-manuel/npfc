"""
Module standardize
===================
This modules is used to standardize molecules and molecular DataFrames.
"""

# standard
import logging
from collections import Counter
from copy import deepcopy
from more_itertools import intersperse
import pkg_resources
# data handling
import json
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdinchi
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize.metal import MetalDisconnector
from rdkit.Chem.MolStandardize.charge import Uncharger
from rdkit.Chem.MolStandardize.normalize import Normalizer
from rdkit.Chem.MolStandardize.tautomer import TautomerCanonicalizer
from rdkit.Chem.Scaffolds import MurckoScaffold
# graph
from networkx import Graph
# docs
from typing import Union
# dev library
from npfc.draw import depict_mol
from npfc import utils
from npfc.filter import Filter


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

DEFAULT_ELEMENTS = {'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'}

# DEFAULT_PROTOCOL = {'tasks': ['filter_empty',
#                               'disconnect_metal',
#                               'clear_mixtures',
#                               'deglycosylate',
#                               'filter_num_heavy_atom',
#                               'filter_molecular_weight',
#                               'filter_num_ring',
#                               'filter_elements',
#                               'clear_isotopes',
#                               'normalize',
#                               'uncharge',
#                               'canonicalize',
#                               'clear_stereo',
#                               ],
#                     'filter_num_heavy_atom': 'num_heavy_atom > 3',
#                     'filter_molecular_weight': 'molecular_weight <= 1000.0',
#                     'filter_num_ring': 'num_ring > 0',
#                     'filter_elements': f'elements in {", ".join(str(x) for x in DEFAULT_ELEMENTS)}',
#                     }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class FullUncharger(Uncharger):
    """A class derived from rdkit.Chem.MolStandardize.charge.Uncharger, so
    instead of attempting to create zwitterions all possible charges are removed
    from the molecule.

    For instance:

    >>> # Uncharger:
    >>> [O-][N+](C)(C)C[O-] -> [O-][N+](C)(C)CO

    >>> # FullUncharger:
    >>> [O-][N+](C)(C)C[O-] -> O[N+](C)(C)CO
    """

    def __init__(self):
        """Create an instance of FullUncharger.

        .. todo:: This will remove charges from -2 to +2 only. This could be improved using more general smarts?
        """
        # some smarts to use to find charges
        self.q_pos_1 = Chem.MolFromSmarts("[*;+]")
        self.q_pos_2 = Chem.MolFromSmarts("[*;+2]")
        self.q_neg_1 = Chem.MolFromSmarts("[*;-]")
        self.q_neg_2 = Chem.MolFromSmarts("[*;-2]")
        logging.debug("Initialized a new FullUncharger object")

    def full_uncharge(self, mol: Mol) -> Mol:
        """Neutralize molecule by adding/removing hydrogens.
        Does not attempt to preserve zwitterions.
        For now takes into account only charges of -2 and +2.

        :param mol: the input molecule
        :return: the uncharged molecule
        """
        logging.debug("Uncharging a molecule")
        mol = deepcopy(mol)
        # Get atom ids for matches
        p = [x[0] for x in mol.GetSubstructMatches(self.q_pos_1)]   # +1
        p += [x[0] for x in mol.GetSubstructMatches(self.q_pos_2)]  # +2
        n = [x[0] for x in mol.GetSubstructMatches(self.q_neg_1)]  # -1
        n += [x[0] for x in mol.GetSubstructMatches(self.q_neg_2)]  # -2
        # remove positive charges
        for atom in [mol.GetAtomWithIdx(x) for x in p]:
            # Remove hydrogen and reduce formal change until neutral or no more hydrogens
            while atom.GetFormalCharge() > 0 and atom.GetNumExplicitHs() > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        # remove negative charges
        for atom in [mol.GetAtomWithIdx(x) for x in n]:
            # Add hydrogen and increase formal change until neutral
            while atom.GetFormalCharge() < 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)

        # clean-up
        mol.ClearComputedProps()
        mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol)

        # mol.ClearComputedProps()  # not testes but might solved the -O-H2 issue
        # mol.UpdatePropertyCache()
        return mol


class Standardizer(Filter):
    """A class for standardizing molecular structures. The standardization itself is based
    on a protocol that the user can modify.

    By default this protocol consists in 15 tasks applied to each molecule invidually:

        1) **initiate_mol**: check if the molecule passed the RDKit conversion
        2) **filter_empty**: filter molecules with empty structures
        3) **disconnect_metal**: break bonds involving metallic atoms, resulting in potentially several molecules per structure.
        4) **clear_mixtures**: retrieve only the "best" molecule from a mixture, which might not always be the largest one.
        5) **deglycosylate**: remove all external sugars-like rings from the molecule and return the remaining non-linear entity.
        6) **filter_num_heavy_atom**: filter molecules with a heavy atom count not in the accepted range. By default: num_heavy_atom > 3.
        7) **filter_molecular_weight**: filter molecules with a molecular weight not in the accepted range. By default: molecular_weight <= 1000.0.
        8) **filter_num_ring**: filter molecules with a number of rings (Smallest Sets of Smallest Rings or SSSR) not in the accepted range. By default: num_ring > 0.
        9) **filter_elements**: filter molecules with elements not considered as medchem. By default: elements in H, B, C, N, O, F, P, S, Cl, Br, I.
        10) **clear_isotopes**: set all atoms to their most common isotope (i.e. 14C becomes 12C which is C).
        11) **normalize**: always write the same functional groups in the same manner.
        12) **uncharge**: remove all charges on a molecule when it is possible. This is different from rdkit.Chem.MolStandardize.charge module as there is no attempt for reaching the zwitterion.
        13) **canonicalize**: enumerate the canonical tautomer.
        14) **clear_stereo**: remove all remaining stereochemistry flags on the molecule.
        15) **reset_mol**: convert forth and back to SMILES format to discard potential residual outdated flags on atoms and bonds.

    Other steps are not part of this protocol but can be executed as well for convenience:

        - **depict**: find the "best" possible 2D depiction of the molecule among Input/rdDepictor/Avalon/CoordGen methods
        - **extract_murcko**: return the Murcko Scaffold from the molecule
        - **clear_side_chains**: remove any exocyclic atom that is not part of a linker
        - **reset_mol**: reset the molecule by converting to and then from smiles

    This results in new columns in the input DataFrame:

        - the 'mol' column: updated structure (only for the protocol)
        - the 'status' column: either passed, filtered or error.
        - the 'task' column: the latest task that was applied to the molecule.

    The standardizer works either on a molecule (method: 'run') or on a DataFrame containing molecules ('run_df').

    In the latter case, the inchikey is computed and can be used for identifying duplicate entries.

    A timeout value is set by default and will be applied to each molecule individually to avoid the process being stuck on marginally difficult cases.
    This value can be set either during the Standardizer object initialization or by defining as an option in the protocol (priority is given to the latter if defined).
    """

    def __init__(self,
                 protocol: str = None,
                 col_mol: str = 'mol',
                 col_id: str = 'idm',
                 elements_medchem: set = DEFAULT_ELEMENTS,
                 timeout: int = 10,
                 ):
        """Create a Standardizer object.

        :param protocol: Either a JSON file or a dictionary. The resultung dictinary needs a 'tasks' key that lists all tasks to be excuted as a list.
        :param col_mol: the column with the molecule for when running the run_df method
        :param col_id: the column with the id for when running the run_df method
        :param filter_duplicates:
        """
        # filter
        super(Standardizer, self).__init__()
        # standardizer
        self._elements_medchem = elements_medchem
        self._col_id = col_id
        self._col_mol = col_mol

        if protocol is None:
            self._protocol = json.load(open(pkg_resources.resource_filename('npfc', 'data/std_mols.json'), 'r'))
        else:
            if isinstance(protocol, str):
                self._protocol = json.load(open(protocol, 'r'))
            else:
                self._protocol = protocol
        # workers
        self.metal_disconnector = MetalDisconnector()
        self.normalizer = Normalizer()
        self.full_uncharger = FullUncharger()
        self.canonicalizer = TautomerCanonicalizer()

        # display information on protocol
        if logging.getLogger().level == logging.DEBUG:
            logging.debug("Successfully instanciated a Standardizer object with protocol:")
            [logging.debug("Task #%s: %s", str(i+1).zfill(2), task) for i, task in enumerate(self._protocol['tasks'])]
            [logging.debug("Option %s %s", opt, value) for opt, value in self._protocol.items() if opt != 'tasks']

    def __repr__(self):
        return f"Standardizer ({len(self._protocol['tasks'])} tasks)"

    def describe(self):
        # init
        pad = max(len(x) for x in self._protocol['tasks'])
        head = 'STANDARIDZER={\n'
        tail = '\n}'
        # define a list of tasks with options in parenthesis
        tasks = list(self._protocol['tasks'])
        for i, task in enumerate(tasks):
            if task in self._protocol.keys():
                opt = self._protocol[task].replace(task.replace('filter_', ''), 'x')
                tasks[i] = f"{task} ({opt})"
        # concatenate all parts and intersperse the tasks with bottow arrows, with step index on the left
        return head + '\n'.join(intersperse('↓'.center(pad+10), [str(i+1).zfill(2).ljust(5) + x.center(pad) for i, x in enumerate(tasks)])) + tail

    @property
    def protocol(self):
        return self._protocol

    @protocol.setter
    def protocol(self, protocol: str):
        # input is a json file => convert it to a dict
        if isinstance(protocol, str):
            utils.check_arg_config_file(protocol)
            with open(protocol) as f:
                protocol = json.load(f)
        # input is a dict
        if isinstance(protocol, dict):
            if 'tasks' not in protocol.keys():
                raise ValueError("invalid protocol format (no 'tasks' key found)")
            elif not isinstance(protocol['tasks'], list) and not isinstance(protocol['tasks'], tuple):
                raise ValueError("invalid protocol format ('tasks' key is neither list or tuple)")
        # update default protocol
        self._protocol.update(protocol)

    @property
    def timeout(self) -> str:
        return self._timeout

    @timeout.setter
    def timeout(self, value: int) -> None:
        if not isinstance(value, int):
            raise ValueError(f"Error! timeout should be a positive int (>1), not '{type(value)}'.")
        elif value < 1:
            raise ValueError(f"Error! timeout should be superior to 1 ({value})")
        self._col_id = value

    @property
    def col_id(self) -> str:
        return self._col_id

    @col_id.setter
    def col_id(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! col_id cannot be '{value}'.")
        self._col_id = value

    @property
    def col_mol(self) -> str:
        return self._col_mol

    @col_mol.setter
    def col_mol(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! col_mol cannot be '{value}'.")
        self._col_mol = value

    @property
    def elements_medchem(self) -> set:
        return self._elements_medchem

    @elements_medchem.setter
    def elements_medchem(self, value: set) -> None:
        if not isinstance(value, set):
            raise ValueError(f"Error! elements_medchem should be a set of strings, not '{value}' ({type(value)}).")
        elif not all([isinstance(v, str) for v in value]):
            raise ValueError(f"Error! elements_medchem should be a set of strings, not '{value}' ({type(value)}).")
        self._elements_medchem = value

    def clear_isotopes(self, mol: Mol) -> Mol:
        """Return a molecule without any isotopes.

        :param mol: the input molecule
        :return: the molecule without isotope
        """
        mol = Mol(mol)
        for a in mol.GetAtoms():
            a.SetIsotope(0)
        return mol

    def clear_mixtures(self, mol: Mol) -> Mol:
        """Return the "best" molecule found in a molecular structure.

        The "best" molecule is determined by the following criteria, sorted by priority:

            1) contains only medchem elements
            2) contains at least one ring
            3) has the largest molecular weight of the mixture

        To summarize:

        .. math::
            medchem > non linear > molecular weight

        So the largest molecule of a mixture might not always be selected, for instance
        a very long aliphatic chain would be dismissed to keep a benzene molecule instead.

        This is implemented in such a way because our fragments used for substructure search contain at least one ring.
        On the contrary, this long aliphatic chain would be kept in a mixture with a non-medchem molecule.

        :param mol: the input molecule(s)
        :return: the best molecule
        """
        submols = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        # no need to look further if we have only one submol!
        if len(submols) < 2:
            return mol
        # otherwise, we have to compare the submols
        # init
        logging.debug("found %s submols", len(submols))
        best_molecular_weight = -1.0  # so we are sure to update this on the first iteration
        best_submol = None
        best_is_medchem = False
        best_is_non_linear = False
        # begin
        for i, submol in enumerate(submols):
            # is_medchem
            is_medchem = self.filter_mol(submol, f'elements in {", ".join(str(x) for x in self.elements_medchem)}')
            is_non_linear = self.filter_mol(submol, "num_rings > 0")
            # molecular_weight
            molecular_weight = Descriptors.ExactMolWt(submol)
            logging.debug("submol #%s: IM=%s, INL=%s, MW=%s", i, is_medchem, is_non_linear, molecular_weight)
            # compare to the current best fragment
            update_best = False
            compute_diff = False  # check which
            # 2 criteria more important than molecular weight: is_medchem > is_non_linear
            if not best_is_medchem and is_medchem:
                update_best = True
            elif best_is_medchem and not is_medchem:
                continue
            elif not best_is_medchem and not is_medchem:
                if not best_is_non_linear and is_non_linear:
                    update_best = True
                elif best_is_non_linear and not is_non_linear:
                    continue
                else:
                    compute_diff = True
            else:  # best_is_medchem and is_medchem
                if not best_is_non_linear and is_non_linear:
                    update_best = True
                elif best_is_non_linear and not is_non_linear:
                    continue
                else:
                    compute_diff = True

            # check molecular_weights only in case of doubt
            if not update_best and compute_diff and molecular_weight > best_molecular_weight:
                update_best = True

            # update best with the properties of the current mol
            if update_best:
                best_is_medchem = is_medchem
                best_is_non_linear = is_non_linear
                best_submol = submol
                best_molecular_weight = molecular_weight

        return best_submol

    def clear_side_chains(self, mol: Mol, debug: bool = False) -> Mol:
        """Clear the side chains of a molecule.

        This method operates in 3 steps:

            1. Remove quickly all atoms in side chains but the one attached to a ring, starting from the terminal atom. (would certainly fail in case of linear molecules)
            2. Iterate over each remaining exocyclic atoms to remove only atoms when it does not break the ring aromaticity. Simple and double bonds can be broken and the atoms in rings which were attached to removed atoms are neutralized.
            3. Remove eventual nitrogen radicals by Smiles editing.

        .. warning:: I found only nitrogen radicals in my dataset, this might be insufficient on a larger scale.

        .. warning:: I found a bug for this molecule 'O=C(O)C1OC(OCC2OC(O)C(O)C(O)C2O)C(O)C(O)C1O', where a methyl remains after processing.

        :param mol: the molecule to simplify
        :return: a simplified copy of the molecule
        """
        # 1st peeling: fast, chunks of terminal chains
        smarts = Chem.MolFromSmarts('[!#1;R0;D1]~[!#1;R0;D{1-2}]')  # terminal exocyclic atom linked to another exocyclic atom, neighbour atom is not allowed more than 2 degrees, so branches (i.e. CC(=O)C) are not cut out
        while mol.HasSubstructMatch(smarts):
            mol = Chem.DeleteSubstructs(mol, smarts)
            mol.ClearComputedProps()
            mol.UpdatePropertyCache()
            Chem.GetSymmSSSR(mol)
        # 2nd peeling: slow, atom per atom of the remaining termninal atoms
        rwmol = Chem.RWMol(mol)
        smarts = Chem.MolFromSmarts('[!#1;R0;D1]')  # remaining terminal exocyclic atoms require cautious handling
        matches = sorted([item for sublist in rwmol.GetSubstructMatches(smarts) for item in sublist], reverse=True)  # reverse order so that remaining atom indices from matches are still valid after removing an atom
        for m in matches:
            try:
                m = m[0]  # should be only single atoms
                rwmol_tmp = deepcopy(rwmol)
                neighbor = rwmol_tmp.GetAtomWithIdx(m).GetNeighbors()[0]  # terminal atom so only 1 neighbor
                rwmol_tmp.RemoveAtom(m)
                neighbor.SetFormalCharge(0)  # neutralize in case of previously quaternary nitrogens
                neighbor.SetNumRadicalElectrons(0)  # remove radicals,this does not work as expected
                Chem.SanitizeMol(rwmol_tmp)  # will fail in case of break in aromaticity
                rwmol = rwmol_tmp  # if it went ok
            except Chem.rdchem.KekulizeException:
                pass  # we should not have tried to remove this atom, so just leave it be

        # I could not figure out how to remove radicals, so I just convert the mol to Smiles and edit the string
        return Chem.MolFromSmiles(Chem.MolToSmiles(rwmol).replace('[N]', 'N').replace('[n]', 'n'))

    def _is_sugar_like(self, ring_aidx: list, mol: Mol):
        """Indicate whether a ring (defined by its atom indices) in a molecule is sugar-like or not.

        Several conditions are to be met for a ring to be considered sugar-like:

            1. size: either 5 or 6 atoms
            2. elements: 1 oxygen and the rest carbons
            3. hybridization: ring atoms need have single bonds only
            4. connection points (next to the ring oxygen): at least 1 has an oxygen as neighbor
            5. subsituents (not next tot the ring oxygen): at least 1/2 (for 5-6-membered rings) have an oxygen as neighbor

        :param ring_aidx: the molecule indices of the ring to investigate
        :param mol: the molecule that contain the ring
        :return: True if the ring complies to the 5 conditions above, False otherwise.
        """
        # ring size: only 5-6 membered rings, rings are already fused when this function is called
        ring_size = len(ring_aidx)
        if ring_size != 5 and ring_size != 6:
            return False

        # access the actual atom objects quickier
        ring_atoms = [mol.GetAtomWithIdx(x) for x in ring_aidx]  # ring atoms are in the same order as ring_aidx

        # atom composition
        elements = [x.GetAtomicNum() for x in ring_atoms]
        element_counter = Counter(elements)
        if not ((ring_size == 5 and element_counter[6] == 4 and element_counter[8] == 1) or (ring_size == 6 and element_counter[6] == 5 and element_counter[8] == 1)):
            return False

        # hybridization of carbon atoms (check if only single bonds attached to the ring)
        carbon_atoms = [x for x in ring_atoms if x.GetAtomicNum() == 6]
        if any([x for x in carbon_atoms if x.GetHybridization() != 4]):  # to check if no H attached in case of the * position
            return False

        # to define connection points and substituents, we first need to identify the position of the ring oxygen
        oxygen_aidx = [x for x in ring_atoms if x not in carbon_atoms][0].GetIdx()  # only 1 oxygen in ring

        # connection points: 1 need at least 1 oxygen as neighbor
        cps = []
        cps_ok = False
        for carbon_atom in carbon_atoms:
            neighbors = carbon_atom.GetNeighbors()
            # if the ring oxygen is next to this atom, this atom is a connection point
            if any([n.GetIdx() == oxygen_aidx for n in neighbors]):
                cps.append(carbon_atom)
                # at least 1 of the connection points has to have an oxygen as side chain
                if any([n.GetAtomicNum() == 8 and n.GetIdx() != oxygen_aidx for n in neighbors]):
                    cps_ok = True
        if not cps_ok:
            return False

        # substituents
        substituents = [x for x in carbon_atoms if x.GetIdx() not in [y.GetIdx() for y in cps]]
        count_oxygens = 0
        for substituent in substituents:
            side_chain_atoms = [x for x in substituent.GetNeighbors() if x.GetIdx() not in ring_aidx]
            if len(side_chain_atoms) > 0:
                if not side_chain_atoms[0].GetAtomicNum() == 8:  # do not check for the degree here because there are connections on substituents too!
                    return False
                count_oxygens += 1
        # at least 1 oxygen for 5-membered rigns and 2 for 6-membered rings
        if (ring_size == 6 and count_oxygens < 2) or (ring_size == 5 and count_oxygens < 1):
            return False

        return True

    def deglycosylate(self, mol: Mol, mode: str = 'run') -> Union[Mol, Graph]:
        """Function to deglycosylate molecules.

        Several rules are applied for removing Sugar-Like Rings (SLRs) from molecules:

            1. Only external SLRs are removed, so a molecule with aglycan-SLR-aglycan is not modified
            2. Only molecules with both aglycans and SLRs are modified (so only SLRs or none are left untouched)
            3. Linear aglycans are considered to be part of linkers and are thus never returned as results
            4. Glycosidic bonds are defined as either O or CO and can be linked to larger linear linker. So from a SLR side, either nothing or only 1 C are allowed before the glycosidic bond oxygen
            5. Linker atoms until the glycosidic bond oxygen atom are appended to the definition of the SLR, so that any extra methyl is also removed.


        .. image:: _images/std_deglyco_algo.svg
            :align: center

        :param mol: the input molecule
        :param mode: either 'run' for actually deglycosylating the molecule or 'graph' for returning a graph of rings instead (useful for presentations or debugging)
        :return: the deglycosylated molecule or a graph of rings
        """

        if len(Chem.GetMolFrags(mol)) > 1:
            raise ValueError("Error! Deglycosylation is designed to work on single molecules, not mixtures!")

        if mode not in ('run', 'graph'):
            raise AttributeError(f"Error! Unauthorized value for parameter 'mode'! ('{mode}')")

        # avoid inplace modifications
        mol = Chem.Mol(mol)

        # define rings
        rings = mol.GetRingInfo().AtomRings()
        rings = utils.fuse_rings(rings)
        # try to deglycosylate only if the molecule has at least 2 rings:
        # - leave linear compounds out
        # - leave sugars in case they are the only ring on the molecule
        if len(rings) < 2:
            return mol

        # annotate sugar-like rings
        are_sugar_like = [self._is_sugar_like(x, mol) for x in rings]
        logging.debug('RINGS: %s', [(rings[i], are_sugar_like[i]) for i in range(len(rings))])
        # remove sugars only when the molecule has some sugar rings and is not entirely composed of sugars
        if not any(are_sugar_like) or all(are_sugar_like):
            return mol
        ring_atoms = set([item for sublist in rings for item in sublist])

        # init sugar graph
        G = Graph()
        # init linkers parts from left and right the glycosidic bond oxygen: one of the side is required to have either C or nothing
        authorized_linker_parts = [[], ['C']]  # R1-OxxxxR2 or R1-COxxxxR2 with xxxx being any sequence of linear atoms (same for R2->R1)

        # define linker atoms as shortest path between 2 rings that do not include other rings
        for i in range(len(rings)):
            ring1 = rings[i]
            for j in range(i+1, len(rings)):
                ring2 = rings[j]
                logging.debug('NEW RING PAIR -- R1: %s; R2: %s', ring1, ring2)

                # shortest path between the two rings that do not include the current rings themselves
                shortest_path = [x for x in Chem.GetShortestPath(mol, ring1[0], ring2[0]) if x not in ring1 + ring2]
                # define the other ring atoms
                other_ring_atoms = ring_atoms.symmetric_difference(set(ring1 + ring2))
                # shortest path for going from the left (ring1) to the right (ring2)
                shortest_path_elements = [mol.GetAtomWithIdx(x).GetSymbol() for x in shortest_path]

                # in case ring1 (left) or/and ring2 (right) is sugar-like, append the side chains left and right
                # to the oxygen to the corresponding ring atoms to avoid left-overs (the O remains is not removed)
                glycosidic_bond = False
                if 'O' in shortest_path_elements:  # not expected to be common enough for a try/catch statement
                    # from the left side
                    aidx_oxygen_left = shortest_path_elements.index('O')  # first O found in list
                    logging.debug('R1 -> R2 -- pos of O: %s; R1 is sugar_like: %s; linker part from R1: %s', aidx_oxygen_left, are_sugar_like[i], shortest_path_elements[:aidx_oxygen_left])
                    if are_sugar_like[i] and shortest_path_elements[:aidx_oxygen_left] in authorized_linker_parts:
                        glycosidic_bond = True
                        ring1 += shortest_path[:aidx_oxygen_left]

                    # from the right side
                    shortest_path_elements.reverse()
                    shortest_path.reverse()
                    aidx_oxygen_right = shortest_path_elements.index('O')  # first O found in list
                    logging.debug('R2 -> R1 -- pos of O: %s; R2 is sugar_like: %s; linker part from R2: %s', aidx_oxygen_right, are_sugar_like[j], shortest_path_elements[:aidx_oxygen_right])
                    if are_sugar_like[j] and shortest_path_elements[:aidx_oxygen_right] in authorized_linker_parts:
                        glycosidic_bond = True
                        ring2 += shortest_path[:aidx_oxygen_right]
                logging.debug('R1 and R2 are linked through a glycosidic bond: %s', glycosidic_bond)

                # in case the 2 rings are directly connected, append a new edge to G
                if not set(shortest_path).intersection(other_ring_atoms):
                    G.add_edge(i, j, atoms=''.join(shortest_path_elements), glycosidic_bond=glycosidic_bond)
                    # annotate nodes with the ring atoms (+ relevent linker atoms) and if they are sugar-like
                    G.nodes[i]['atoms'] = ring1
                    G.nodes[i]['sugar_like'] = are_sugar_like[i]
                    G.nodes[j]['atoms'] = ring2
                    G.nodes[j]['sugar_like'] = are_sugar_like[j]

        # draw the graph
        if mode == 'graph':
            # colormap_nodes = [(0.7,0.7,0.7) if x['sugar_like'] else (1,0,0) for i, x in G.nodes(data=True)]
            # return draw.fc_graph(G, colormap_nodes=colormap_nodes)
            return G

        # iterative recording of terminal sugar rings (atoms) that are linked with a glycosidic bond
        ring_atoms_to_remove = []
        nodes_to_remove = [node for node in G.nodes(data=True) if node[1]['sugar_like'] and G.degree(node[0]) == 1 and list(G.edges(node[0], data=True))[0][2]['glycosidic_bond']]
        while len(nodes_to_remove) > 0:
            # record atoms indices to remove from the molecule
            [ring_atoms_to_remove.append(n[1]['atoms']) for n in nodes_to_remove]
            # remove nodes from current layer for next iteration
            [G.remove_node(n[0]) for n in nodes_to_remove]
            nodes_to_remove = [node for node in G.nodes(data=True) if node[1]['sugar_like'] and G.degree(node[0]) == 1 and list(G.edges(node[0], data=True))[0][2]['glycosidic_bond']]
        logging.debug('Ring atoms to remove: %s', ring_atoms_to_remove)

        # edit the molecule
        if ring_atoms_to_remove:
            # flatten the atom indices of each ring to remove in reverse order so that atom indices do not change when removing atoms
            ring_atoms_to_remove = sorted([item for sublist in ring_atoms_to_remove for item in sublist], reverse=True)
            emol = Chem.EditableMol(mol)
            [emol.RemoveAtom(x) for x in ring_atoms_to_remove]
            mol = emol.GetMol()
            logging.debug('Obtained fragments: %s', Chem.MolToSmiles(mol))

        # clean-up
        frags = Chem.GetMolFrags(mol, asMols=True)
        # avoid counting the number of rings in each fragment if only 1 fragment left anyway
        if len(frags) == 1:
            logging.debug('Only one fragment obtained, returning it')
            return frags[0]
        # the substituents of the deleted terminal sugar-like rings remain in the structure,
        # these are obligatory linear because they were not in the graph,
        # so we just have to retrieve the one fragment that is not linear
        logging.debug('Returning only the non-linear obtained fragment')
        return [x for x in frags if Descriptors.rdMolDescriptors.CalcNumRings(x) > 0][0]

    def _run(self, mol: Mol) -> tuple:
        """Helper function for run.
        Contains all tasks defined within the protocol. Since some operations are
        expansive and could last a very long time for complex molecules (normalize, canonicalize),
        a timeout value is set globally. The run function is the one that can catch the exception
        raised by timeouts.

        :param mol: the input molecule
        :return: a tuple containing the molecule, its status and the further task name it reached
        """

        # initiate_mol
        if mol is None:
            return (mol, 'error', 'initiate_mol')

        # begin protocol
        mol = deepcopy(mol)  # do not modify the molecule in place
        for task in self._protocol['tasks']:
            # filter_empty
            if task == 'filter_empty':
                try:
                    if not mol.GetNumAtoms():
                        return ('', 'filtered', task)
                except ValueError:
                    return ('', 'error', task)

            # disconnect_metal
            elif task == 'disconnect_metal':
                try:
                    mol = self.metal_disconnector.disconnect(mol)
                except ValueError:
                    return (mol, 'error', 'disconnect_metal')

            # clear_mixtures
            elif task == 'clear_mixtures':
                try:
                    mol = self.clear_mixtures(mol)
                except ValueError:
                    return (mol, 'error', task)

            # deglycosylate
            elif task == 'deglycosylate':
                try:
                    mol = self.deglycosylate(mol)
                except ValueError:
                    return (mol, 'error', task)

            # filters
            elif task.startswith('filter_'):  # filter empty is tried before
                try:
                    if not self.filter_mol(mol, self._protocol[task]):
                        return (mol, 'filtered', task)
                except ValueError:
                    return (mol, 'error', task)

            # sanitize
            elif task == 'sanitize':
                try:
                    Chem.SanitizeMol(mol)
                except ValueError:
                    return (mol, 'error', task)

            # clear_isotopes
            elif task == 'clear_isotopes':
                try:
                    mol = self.clear_isotopes(mol)
                except ValueError:
                    return (mol, 'error', task)

            # normalize
            elif task == 'normalize':
                try:
                    mol = self.normalizer.normalize(mol)
                except ValueError:
                    return (mol, 'error', task)

            # uncharge
            elif task == 'uncharge':
                try:
                    mol = self.full_uncharger.full_uncharge(mol)
                except ValueError:
                    return (mol, 'error', task)

            # canonicalize
            elif task == 'canonicalize':
                # canonicalize
                try:
                    mol = self.canonicalizer.canonicalize(mol)
                except (ValueError, RuntimeError):
                    return (mol, 'error', task)

            # clear_stereo
            elif task == 'clear_stereo':
                try:
                    rdmolops.RemoveStereochemistry(mol)
                except ValueError:
                    return (mol, 'error', task)

            # extract Murcko Scaffolds
            elif task == 'extract_murcko':
                try:
                    mol = MurckoScaffold.GetScaffoldForMol(mol)
                except ValueError:
                    return (mol, 'error', task)

            # clear side chains
            elif task == 'clear_side_chains':
                try:
                    mol = self.clear_side_chains(mol)
                except ValueError:
                    return (mol, 'error', task)

            elif task == 'depict':
                try:
                    mol = depict_mol(mol)
                except ValueError:
                    return (mol, 'error', task)

            elif task == 'reset_mol':
                try:
                    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
                except ValueError:
                    return (mol, 'error', task)

            # something else?
            else:
                raise ValueError(f"Unknown task: {task}")

        # a molecule that passed all the protocole!
        return (mol, 'passed', 'standardize')

    def run(self, mol: Mol, timeout: int = 10) -> tuple:
        """Execute the standardization protocol on a molecule.
        Molecule that exceed the timeout value are filtered with a task='timeout'.

        As a final step of the protocol, InChiKeys ('inchikey') are computed for identifying molecules.

        :param mol: the input molecule
        :param timeout: the maximum number of seconds for processing a molecule
        :return: a tuple containing the molecule, its status and the further task name it reached
        """
        with utils.timeout(timeout):
            return self._run(mol)

        # in case of timeout
        return (mol, 'filtered', 'timeout')


    def run_df(self, df: DataFrame) -> tuple:
        """Apply the standardization protocol on a DataFrame, with the possibility of directly filtering duplicate entries as well.
        This can be very useful as the standardization process can expose duplicate entries due to salts removal, neutralization,
        canonical tautomer enumeration, and stereochemistry centers unlabelling

        If a reference file is specified, duplicate removals becomes possible accross chunks.


        :param df: the input DataFrame
        :param timeout: the maximum number of seconds for processing a molecule
        :return: three DataFrames separated by status:

            - passed
            - filtered
            - error

        .. note:: As a side effect, the output DataFrames get indexed by idm. The 'inchikey' col is not returned, but the values can be accessed using the reference file.

        :param df: The DataFrame with molecules to standardize
        :param return: a tuple of 3 DataFrames: standardized, filtered and error.
        """
        # run standardization protocol
        df.index = df[self.col_id]
        df.loc[:, self.col_mol], df.loc[:, 'status'], df.loc[:, 'task'] = zip(*df[self.col_mol].map(self.run))
        # flag eventual None molecules at the end of the pipeline for filtering out
        df['status'] = df.apply(lambda x: x['status'] if x['mol'] is not None else 'error', axis=1)
        df['task'] = df.apply(lambda x: x['task'] if x['mol'] is not None else 'filter_empty_final', axis=1)
        df['mol'] = df['mol'].map(lambda x: x if x is not None else '')
        # do not apply filter duplicates on molecules with errors or that were already filtered for x reasons
        df_error = df[df['status'] == 'error']
        df_filtered = df[df['status'] == 'filtered']
        df = df[df['status'].str.contains('passed')]
        df = df.copy()  # only way I found to suppress pandas warnings in a "clean" way
        # compute InChiKeys
        df.loc[:, 'inchikey'] = df.loc[:, self.col_mol].map(rdinchi.MolToInchiKey)

        # tuple of dataframes
        return (df, df_filtered, df_error)
