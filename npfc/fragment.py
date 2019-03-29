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

    def classify_fragment_combination(self,
                                      mol: Mol,
                                      aidxf1: set,
                                      aidxf2: set,
                                      cutoff: int = 3) -> dict:
        """Classify a fragment combination found in a molecule as a dictionary
        with category, type and subtype values.

        Possible classifications are:

        - fusion
            - spiro fsp)
            - edge (fed)
            - bridged (fbr)
            - unkown (fun)
            - false_positive
                - substructure (ffs)
                - overlap (ffo)
        - connection
            - monopodal (cmo)
            - bipodal
                - spiro (cbs)
                - edge (cbe)
                - bridged (cbb)
            - tripodal
                - spiro (cts)
                - edge (cte)
                - bridged (ctb)
            - unknown
                - spiro (cus)
                - edge (cue)
                - bridged (cub)
            - false_positive
                - cutoff (cfc)

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
        if aidxf1.issubset(aidxf2):
            return {'category': 'fusion', 'type': 'false_positive', 'subtype': 'substructure', 'abbrev': 'ffs'}
        elif aidxf2.issubset(aidxf1):
            return {'category': 'fusion', 'type': 'false_positive', 'subtype': 'substructure', 'abbrev': 'ffs'}
        if len(aidx_fused) > 0:
            category = 'fusion'
            if len(aidx_fused) == 1:
                return {'category': category, 'type': 'spiro', 'subtype': '', 'abbrev': 'fsp'}
            elif len(aidx_fused) == 2:
                return {'category': category, 'type': 'edge', 'subtype': '', 'abbrev': 'fed'}
            elif len(aidx_fused) == 3:
                return {'category': category, 'type': 'bridged', 'subtype': '', 'abbrev': 'fbr'}
            else:
                sssr = mol.GetRingInfo().AtomRings()  # smallest sets of smallest rings
                for aidxr in sssr:
                    # if at least one ring is completely present in the overlap between fragments,
                    # then it's a false positive due to the fragments overlap qnd not a combination.
                    if set(aidxr).issubset(aidx_fused):
                        return {'category': category, 'type': 'false_positive', 'subtype': 'overlap', 'abbrev': 'ffo'}
                # something unknown with > 3 fused atoms!
                return {'category': category, 'type': 'unknown', 'subtype': '', 'abbrev': 'fun'}
        else:
            category = 'connection'
            logging.debug(f"category: {category}")
            shortest_path_between_frags = self.get_shortest_path_between_frags(mol, aidxf1, aidxf2)
            logging.debug(f"shortest_path_between_frags: {shortest_path_between_frags} ({len(shortest_path_between_frags)})")
            if len(shortest_path_between_frags) - 2 > cutoff:  # begin and end atoms are in the shortest path but should not be considered for cutoff
                logging.debug(f"shortest_path_between_frags is greater than cutoff ({len(shortest_path_between_frags) - 2} > {cutoff})")
                return {'category': category, 'type': 'false_positive', 'subtype': 'cutoff', 'abbrev': 'cfc'}
            intermediary_rings = self.get_rings_between_two_fragments(mol, aidxf1, aidxf2)
            if len(intermediary_rings) == 0:
                return {'category': category, 'type': 'monopodal', 'subtype': '', 'abbrev': 'cmo'}
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
        :return: the dictionary specifying fragment combination category, type, subtype and abbrev
        """
        for ir in intermediary_rings:
            intersect_1 = ir.intersection(aidxf1)
            intersect_2 = ir.intersection(aidxf2)
            logging.debug(f"ir:{ir} => intersect_1: {intersect_1} ({len(intersect_1)}), intersect_2: {intersect_2} ({len(intersect_2)})")
            if len(intersect_1) == 1 or len(intersect_2) == 1:
                return {'category': category, 'type': type, 'subtype': 'spiro', 'abbrev': category[0] + type[0] + 's'}
            elif len(intersect_1) == 3 or len(intersect_2) == 3:
                return {'category': category, 'type': type, 'subtype': 'bridged', 'abbrev': category[0] + type[0] + 'b'}
            elif len(intersect_1) > 3 or len(intersect_2) > 3:
                return {'category': category, 'type': type, 'subtype': 'unknown', 'abbrev': category[0] + type[0] + 'u'}
        # edge otherwise
        return {'category': category, 'type': type, 'subtype': 'edge', 'abbrev': category[0] + type[0] + 'e'}

    def classify_fragment_combinations(self,
                                       df_mols: DataFrame,
                                       df_aidxf: DataFrame,
                                       cutoff: int = 3) -> DataFrame:
        """Return a DataFrame with all fragment combination categories for a given set of
        molecules and fragment atom indices obtained by substructure search.
        For more details about category, type and subtype, see doc in method classify_fragment_combination.

        The output DataFrame contains 8 columns decribing each fragment combination:

        1) idm: the id of the molecule
        2) idf1: the id of fragment 1
        3) idf2: the id of fragment 2
        4) abbrev: a 3-letter code indicating category, type and subtype
        5) category
        6) type
        7) subtype
        8) aidxf1: the atom indices of fragment 1 found in the molecule
        9) aidxf2: the atom indices of fragment 2 found in the molecule

        Fragments with a number of intermediary atoms higher than defined cutoff
        are labelled as false positives.

        :param df_mols: the input DataFrame with molecules
        :param df_aidxf: the input DataFrame with substructure matches
        :param cutoff: the maximum number of intermediary atoms between 2 fragments
        :return: a DataFrame with all fragment combination classifications
        """
        ds_fcc = []
        if 'idm' in df_mols.columns:
            logging.debug("Reindexing df with 'idm'")
            df_mols.index = df_mols['idm']
        # labelling idxf
        df_aidxf['aidxf_str'] = df_aidxf['aidxf'].map(str)  # sets are an unhashable type...
        df_aidxf['idxf'] = df_aidxf.groupby(['idm', 'idf', 'aidxf_str']).grouper.group_info[0]

        # classify fragment combinations
        for gid, g in df_aidxf.groupby('idm'):
            mol = df_mols.loc[gid]['mol']
            # mol = df_mols[df_mols['idm'] == gid]['mol'].iloc[0]
            for i in range(len(g)):
                row_f1 = g.iloc[i]
                aidxf1 = row_f1['aidxf']
                idf1 = row_f1['idf']
                idxf1 = row_f1['idxf']
                for j in range(i+1, len(g)):
                    row_f2 = g.iloc[j]
                    aidxf2 = row_f2['aidxf']
                    idf2 = row_f2['idf']
                    idxf2 = row_f2['idxf']
                    logging.debug(f"Classifying m={gid}, f1={idf1}, f2={idf2}")
                    d_fcc = self.classify_fragment_combination(mol, aidxf1, aidxf2, cutoff=cutoff)

                    # record fragment combination
                    d_fcc['idm'] = gid
                    d_fcc['idf1'] = idf1
                    d_fcc['idxf1'] = idxf1
                    d_fcc['fid1'] = str(idf1) + ":" + str(idxf1)
                    d_fcc['idf2'] = idf2
                    d_fcc['idxf2'] = idxf2
                    d_fcc['fid2'] = str(idf2) + ":" + str(idxf2)
                    d_fcc['aidxf1'] = aidxf1
                    d_fcc['aidxf2'] = aidxf2
                    ds_fcc.append(d_fcc)
        # dataframe with columns in given order
        return DataFrame(ds_fcc, columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', 'aidxf1', 'aidxf2', 'fid1', 'fid2'])

    def clean(self, df_fcc):
        """Clean a df_fcc by removing false positives such as substructures and
        cutoff combinations.
        Currently, the cluster terminated the vast majority of my jobs because of
        memory usage, so I'm trying to split up the different tasks in different
        scripts.
        """

        # clean the data

        logging.debug("Now cleaning fragment combinations")

        # drop cutoff combinations
        logging.debug(f"Removing cutoff connections from fragment combinations")
        df_fcc = df_fcc[df_fcc['abbrev'] != 'cfc']

        # drop fragments combinations paired with a substructure
        df_substructures = df_fcc[df_fcc['abbrev'] == 'ffs']  # all the substructures in the whole dataframe
        logging.debug(f"Number of substructures found in df_fcc: {len(df_substructures.index)}/{len(df_fcc.index)}")
        if len(df_substructures) > 0:
            logging.debug(f"Removing substructures from fragment combinations")
            for gid, g in df_fcc[df_fcc['idm'].isin(df_substructures['idm'])].groupby('idm'):  # iterate only on the groups with at least one substructure
                idf_to_remove = []
                for row in g[g['abbrev'] == 'ffs'].itertuples():
                    if len(row[10]) > len(row[11]):
                        idf_to_remove.append(row[5])
                    else:
                        idf_to_remove.append(row[3])
            df_fcc = df_fcc[(~df_fcc['idxf1'].isin(idf_to_remove)) & (~df_fcc['idxf2'].isin(idf_to_remove))]
        if len(df_fcc.index) == 0:
            logging.debug("No fragment remaining for mapping!")
            return None
        logging.debug(f"Remaining number of fragment combinations: {len(df_fcc.index)}")

        return df_fcc

    def map_frags(self, df_fcc: DataFrame, min_frags: int = 2, max_frags: int = 5, max_overlaps: int = 5) -> DataFrame:
        """
        This method process a fragment combinations computed with classify_fragment_combinations
        and return a new DataFrame with a fragment map for each molecule.

        This fragment map is a single line string representation of the fragment connectivity
        within a molecule and follows following syntax:

            >>> f1[abbrev1]f2-f1[abbrev2]f3-f2[abbrev3]f3

        No applying any limit of the max number of frags might have been what caused
        crashed due to memory usage on the cluster.
        Nope, this still happens now.
        The real reason was because of very high numbers of overlapping combinations
        in some molecules (366 so 2^366 graphs!)
        """
        # split by overlaps

        logging.debug(f"Mapping fragments")

        ds_map = []
        for gid, g in df_fcc.groupby('idm'):
            # entries with an overlap
            overlaps = g[g['abbrev'] == 'ffo']
            noverlaps = len(overlaps.index)
            logging.debut(f"Number of overlaps found for molecule {gid}: {noverlaps}")
            if noverlaps > max_overlaps:
                logging.debug(f"Too many overlap combinations ({noverlaps}), discarding molecule '{gid}'")
                continue
            if len(overlaps.index) > 0:  # the code below could certainly be improved, but this case should not happen too often
                # remove these from the current group
                g = g[g['abbrev'] != 'ffo']
                # get fragment ids of the invariant parts of alternative paths
                common = g[(~g['fid1'].isin(overlaps['fid1'])) & (~g['fid2'].isin(overlaps['fid2']))]
                common_combinations = set()
                for rowid, row in common.iterrows():
                    common_combinations.add(row['fid1'])
                    common_combinations.add(row['fid2'])
                common_combinations = list(common_combinations)
                # get the fragment ids of the variant parts of alternative paths
                alt_combinations = []
                for rowid, row in overlaps.iterrows():
                    alt_combinations.append([row['fid1'], row['fid2']])
                # get all possible paths
                alt_combinations = [list(x) + common_combinations for x in list(product(*alt_combinations))]

                dfs_fcc_clean = []
                for alt in alt_combinations:
                    df_alt = g[(g['fid1'].isin(alt)) | (g['fid2'].isin(alt))]
                    dfs_fcc_clean.append(df_alt)
            else:
                dfs_fcc_clean = [g]

            # fragment map string representation
            for df_fcc_clean in dfs_fcc_clean:
                frag_map_str = '-'.join(list(df_fcc_clean['fid1'].map(str) + "[" + df_fcc_clean['abbrev'] + "]" + df_fcc_clean['fid2'].map(str)))
                frags = list(df_fcc_clean['fid1'].map(str).values) + list(df_fcc_clean['fid2'].map(str).values)
                nfrags = len(frags)
                frags_u = list(set(frags))
                nfrags_u = len(frags_u)
                if nfrags_u < min_frags:
                    logging.debug(f"Too few unique fragment occurrences, discarding graph of n={nfrags_u} for molecule: '{gid}'")
                    continue
                elif nfrags_u > max_frags:
                    logging.debug(f"Too many unique fragment occurrences, discarding graph of n={nfrags_u} for molecule: '{gid}'")
                    continue
                comb = list(df_fcc_clean['abbrev'].values)
                ncomb = len(comb)
                comb_u = list(set(comb))
                ncomb_u = len(comb_u)
                ds_map.append({'idm': gid, 'map_str': frag_map_str, 'nfrags': nfrags, 'nfrags_u': nfrags_u, 'ncomb': ncomb, 'ncomb_u': ncomb_u, 'frags': frags, 'frags_u': frags_u, 'comb': comb, 'comb_u': comb_u})

        return DataFrame(ds_map, columns=['idm', 'map_str', 'nfrags', 'nfrags_u', 'ncomb', 'ncomb_u', 'frags', 'frags_u', 'comb', 'comb_u'])
