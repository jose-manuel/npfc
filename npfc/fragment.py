"""
Module fragment
===============
This modules contains two classes:

    - Matcher for substructure search
    - Classifier for classifying fragment combinations
"""

# standard
import logging
from itertools import product
# data science
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import AllChem

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Matcher:
    """Create a Matcher object."""

    def __init__(self):
        pass

    def run(self, df_mols: DataFrame, df_frags: DataFrame,
            col_mol_mols: str = 'mol', col_mol_frags: str = 'mol') -> DataFrame:
        """Create a DataFrame recording every substructure (fragment) match in the
        input molecule DataFrame.

        A record is composed of 4 fields:

        1) idm: the id of the molecule (rowid from df_mols)
        2) idf: the id of the fragment (rowid from df_frags)
        3) aidxf: the atom indices of the fragment found in the molecule
        4) mol_perc: the percentage of the molecule the fragment represents (based on hac)

        :param df_mols: the input DataFrame with molecules
        :param df_frags: the input DataFrame with fragments to use for substructure search
        :param col_mol_frags: the input DataFrame column name with the fragments
        :return: the substructure matches as a DataFrame

        .. note:: Rowids labelled idm are used for recording of combinations: mols => idm, frags => idf
        """
        # init
        d = {}
        d['idm'] = []
        d['idf'] = []
        d['aidxf'] = []
        d['mol_perc'] = []  # proportion of the molecule the substructure represents
        # begin
        for idm, rowm in df_mols.iterrows():
            hac = rowm[col_mol_mols].GetNumAtoms()
            for idq, rowq in df_frags.iterrows():
                matches = rowm[col_mol_mols].GetSubstructMatches(rowq[col_mol_frags])
                for m in matches:
                    d['idm'].append(idm)
                    d['idf'].append(idq)
                    d['aidxf'].append(set(m))  # set for intersections later
                    d['mol_perc'].append(round(len(m)/hac, 2) * 100)
        return DataFrame(d)


class CombinationClassifier:
    """A class useful for classifying fragment combinations."""

    def __init__(self):
        """Create an instance of FragmentCombinationClassifier."""
        pass

    def get_rings_between_two_fragments(self, mol: Mol, aidxf1: set, aidxf2: set) -> list:
        """Returns the atom indices of every ring that connects two fragments together, defined by atom indices.

        :param mol: the input molecule
        :param aidxf1: the atom indices of the first fragment found in the molecule
        :param aidxf2: the atom indices of the second fragment found in the molecule
        :return: a list of intermediary rings between both fragments and defined by atom indices
        """
        ri = mol.GetRingInfo()
        intermediary_rings = []
        for i, ar in enumerate(ri.AtomRings()):
            sar = set(ar)
            if sar.intersection(aidxf1) and sar.intersection(aidxf2):
                intermediary_rings.append(sar)

        return intermediary_rings

    def get_shortest_path_between_frags(self, mol: Mol, aidxf1: set, aidxf2: set) -> tuple:
        """Return the shortest path within a molecule between two fragments defined by atom indices.
        First and last atom indices are part of respectively fragment 1 and fragment 2, so they should not
        be considered when estimating the distance between fragments.

        (i.e. distance = len(shortest_path) - 2)

        :param mol: The input molecule.
        :param aidxf1: the atom indices of the first fragment found in the molecule
        :param aidxf2: the atom indices of the second fragment found in the molecule
        :return: the atom indices of the shortest path between both fragments. The first index is the attachment point from fragment 1 whereas the last index is the attachment point from fragment 2
        """
        # 1/ compute every pairwise atom combination between both fragments
        pairwise_combinations = product(tuple(aidxf1), tuple(aidxf2))
        logging.debug(f"pairwise_combinations: {pairwise_combinations}")
        # 2/ for each of those, compute the shortest path possible
        all_paths = [AllChem.GetShortestPath(mol, pc[0], pc[1]) for pc in pairwise_combinations]
        logging.debug(f"all_paths: {all_paths}")
        # 3/ return one of the shortest pathes
        return min(all_paths, key=lambda x: len(x))

    def classify_fragment_combination(self, mol: Mol, aidxf1: set, aidxf2: set, cutoff: int = 3) -> dict:
        """Classify a fragment combination found in a molecule as a dictionary
        with category, type and subtype values.

        Possible classifications are:

        - fusion
            - spiro
            - edge
            - bridged
            - unkown
            - false_positive
                - substructure
        - connection
            - monopodal
            - bipodal
                - spiro
                - edge
                - bridged
            - tripodal
                - spiro
                - edge
                - bridged
            - unknown
                - spiro
                - edge (##TODO: find an example molecule for this!)
                - bridged
            - false_positive
                - cutoff

        :param mol: the input molecule
        :param aidxf1: the atom indices of the first fragment found in the molecule
        :param aidxf2: the atom indices of the second fragment found in the molecule
        :return: the dictionary specifying fragment combination category, type and subtype
        """
        logging.debug(f"cuttoff: {cutoff}")
        # get atoms in common in fragment 1 and fragment 2
        aidx_fused = aidxf1.intersection(aidxf2)
        logging.debug(f"aidx_fused: {aidx_fused}")
        # in case of overlapping fragments in the queries, we get overlapping matches
        if aidxf1.issubset(aidxf2) or aidxf2.issubset(aidxf1):
            return {'category': 'fusion', 'type': 'false_positive', 'subtype': 'substructure'}
        if len(aidx_fused) > 0:
            category = 'fusion'
            if len(aidx_fused) == 1:
                return {'category': category, 'type': 'spiro', 'subtype': ''}
            elif len(aidx_fused) == 2:
                return {'category': category, 'type': 'edge', 'subtype': ''}
            elif len(aidx_fused) == 3:
                return {'category': category, 'type': 'bridged', 'subtype': ''}
            else:
                sssr = mol.GetRingInfo().AtomRings()  # smallest sets of smallest rings
                for aidxr in sssr:
                    # if at least one ring is completely present in the overlap between fragments,
                    # then it's a false positive due to the fragments overlap qnd not a combination.
                    if set(aidxr).issubset(aidx_fused):
                        return {'category': category, 'type': 'false_positive', 'subtype': 'overlap'}
                # something unknown with > 3 fused atoms!
                return {'category': category, 'type': 'unknown', 'subtype': ''}
        else:
            category = 'connection'
            logging.debug(f"category: {category}")
            shortest_path_between_frags = self.get_shortest_path_between_frags(mol, aidxf1, aidxf2)
            logging.debug(f"shortest_path_between_frags: {shortest_path_between_frags} ({len(shortest_path_between_frags)})")
            if len(shortest_path_between_frags) - 2 > cutoff:  # begin and end atoms are in the shortest path but should not be considered for cutoff
                logging.debug(f"shortest_path_between_frags is greater than cutoff ({len(shortest_path_between_frags) - 2} > {cutoff})")
                return {'category': category, 'type': 'false_positive', 'subtype': 'cutoff'}
            intermediary_rings = self.get_rings_between_two_fragments(mol, aidxf1, aidxf2)
            if len(intermediary_rings) == 0:
                return {'category': category, 'type': 'monopodal', 'subtype': ''}
            else:
                intermediary_rings = self.get_rings_between_two_fragments(mol, aidxf1, aidxf2)
                if len(intermediary_rings) == 1:
                    type = 'bipodal'  # 1 intermediary ring, which are defined by intersection of aidx, so at least always 1 for each fragment!
                    return self._get_fcc_subtype(category, type, aidxf1, aidxf2, intermediary_rings)
                elif len(intermediary_rings) == 2:
                    type = 'tripodal'
                    return self._get_fcc_subtype(category, type, aidxf1, aidxf2, intermediary_rings)
                else:
                    type = 'unknown'
                    return self._get_fcc_subtype(category, type, aidxf1, aidxf2, intermediary_rings)

    def _get_fcc_subtype(self, category: str, type: str, aidxf1: set, aidxf2: set, intermediary_rings: list) -> tuple:
        """Return the subtype (spiro, edge, bridged) for bipodal, tripodal and unknown connections.

        :param category: the fragment combination category (connection)
        :param type: the fragment combination type (bipodal, tripodal or unknown)
        :param aidxf1: the atom indices of the first fragment found in the molecule
        :param aidxf2: the atom indices of the second fragment found in the molecule
        :param intermediary_rings: the list of intermediary rings between both fragments and defined by atom indices
        :return: the dictionary specifying fragment combination category, type and subtype
        """
        for ir in intermediary_rings:
            intersect_1 = ir.intersection(aidxf1)
            intersect_2 = ir.intersection(aidxf2)
            logging.debug(f"ir:{ir} => intersect_1: {intersect_1} ({len(intersect_1)}), intersect_2: {intersect_2} ({len(intersect_2)})")
            if len(intersect_1) == 1 or len(intersect_2) == 1:
                return {'category': category, 'type': type, 'subtype': 'spiro'}
            elif len(intersect_1) == 3 or len(intersect_2) == 3:
                return {'category': category, 'type': type, 'subtype': 'bridged'}
            elif len(intersect_1) > 3 or len(intersect_2) > 3:
                return {'category': category, 'type': type, 'subtype': 'unknown'}
        # edge otherwise
        return {'category': category, 'type': type, 'subtype': 'edge'}

    def classify_fragment_combinations(self, df_mols: DataFrame, df_aidxf: DataFrame) -> DataFrame:
        """Return a DataFrame with all fragment combination categories for a given set of
        molecules and fragment atom indices obtained by substructure search.
        For more details about category, type and subtype, see doc in method classify_fragment_combination.

        The output DataFrame contains 8 columns decribing each fragment combination:

        1) idm: the id of the molecule
        2) idf1: the id of fragment 1
        3) idf2: the id of fragment 2
        4) category
        5) type
        6) subtype
        7) aidxf1: the atom indices of fragment 1 found in the molecule
        8) aidxf2: the atom indices of fragment 2 found in the molecule

        :param df_mols: the input DataFrame with molecules
        :param df_aidxf: the input DataFrame with substructure matches
        :return: a DataFrame with all fragment combination classifications
        """
        ds_fcc = []
        if 'idm' in df_mols.columns:
            logging.debug("Reindexing df with 'idm'")
            df_mols.index = df_mols['idm']
        for gid, g in df_aidxf.groupby('idm'):
            mol = df_mols.loc[gid]['mol']
            # mol = df_mols[df_mols['idm'] == gid]['mol'].iloc[0]
            for i in range(len(g)):
                aidxf1 = g.iloc[i]['aidxf']
                idf1 = g.iloc[i]['idf']
                for j in range(i+1, len(g)):
                    aidxf2 = g.iloc[j]['aidxf']
                    logging.debug(f"Classifying m={gid}, f1={idf1}, f2={g.iloc[j]['idf']}")
                    d_fcc = self.classify_fragment_combination(mol, aidxf1, aidxf2)

                    # record fragment combination
                    d_fcc['idm'] = gid
                    d_fcc['idf1'] = idf1
                    d_fcc['idf2'] = g.iloc[j]['idf']
                    d_fcc['aidxf1'] = aidxf1
                    d_fcc['aidxf2'] = aidxf2
                    ds_fcc.append(d_fcc)
        # dataframe with columns in given order
        return DataFrame(ds_fcc, columns=['idm', 'idf1', 'idf2', 'category', 'type', 'subtype', 'aidxf1', 'aidxf2'])
