"""
Module fragment
===============
This modules contains two classes:

    - Matcher for substructure search
    - Classifier for classifying fragment combinations
"""

# standard
import logging
import itertools
import json
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
            for idf, rowq in df_frags.iterrows():
                matches = rowm[col_mol_mols].GetSubstructMatches(rowq[col_mol_frags])
                for m in matches:
                    d['idm'].append(idm)
                    d['idf'].append(idf)
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
        pairwise_combinations = itertools.product(tuple(aidxf1), tuple(aidxf2))
        # 2/ for each of those, compute the shortest path possible
        all_paths = [AllChem.GetShortestPath(mol, pc[0], pc[1]) for pc in pairwise_combinations]
        # logging.debug(f"Looking for the shortest path shortest path among these:")
        # [logging.debug(f"Path ({str(i).zfill(3)}): {p}") for i, p in enumerate(all_paths)]
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
            - other (fot)
            - false_positive
                - substructure (ffs)
                - overlap (ffo)
        - connection
            - monopodal (cmo)
            - bipodal
                - spiro (cbs)
                - edge (cbe)
                - bridged (cbb)
                - other (cbo)
            - tripodal
                - spiro (cts)
                - edge (cte)
                - bridged (ctb)
                - other (cto)
            - other
                - spiro (cos)
                - edge (coe)
                - bridged (cob)
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
            abbrev = 'ffs'
            logging.debug(f"Classification: {abbrev}")
            return {'category': 'fusion', 'type': 'false_positive', 'subtype': 'substructure', 'abbrev': abbrev}
        elif aidxf2.issubset(aidxf1):
            abbrev = 'ffs'
            logging.debug(f"Classification: {abbrev}")
            return {'category': 'fusion', 'type': 'false_positive', 'subtype': 'substructure', 'abbrev': abbrev}
        if len(aidx_fused) > 0:
            category = 'fusion'
            if len(aidx_fused) == 1:
                abbrev = 'fsp'
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': 'spiro', 'subtype': '', 'abbrev': abbrev}
            elif len(aidx_fused) == 2:
                abbrev = 'fed'
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': 'edge', 'subtype': '', 'abbrev': abbrev}
            elif 3 <= len(aidx_fused) <= 5:
                abbrev = 'fbr'
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': 'bridged', 'subtype': '', 'abbrev': 'fbr'}
            else:
                sssr = mol.GetRingInfo().AtomRings()  # smallest sets of smallest rings
                for aidxr in sssr:
                    # if at least one ring is completely present in the overlap between fragments,
                    # then it's a false positive due to the fragments overlap qnd not a combination.
                    if set(aidxr).issubset(aidx_fused):
                        abbrev = 'ffo'
                        logging.debug(f"Classification: {abbrev}")
                        return {'category': category, 'type': 'false_positive', 'subtype': 'overlap', 'abbrev': abbrev}
                # something unknown with > 5 fused atoms!
                return {'category': category, 'type': 'unknown', 'subtype': '', 'abbrev': 'fot'}
        else:
            category = 'connection'
            logging.debug(f"category: {category}")
            shortest_path_between_frags = self.get_shortest_path_between_frags(mol, aidxf1, aidxf2)
            logging.debug(f"shortest_path_between_frags: n={len(shortest_path_between_frags)-2} {shortest_path_between_frags}")
            if len(shortest_path_between_frags) - 2 > cutoff:  # begin and end atoms are in the shortest path but should not be considered for cutoff
                abbrev = 'cfc'
                logging.debug(f"Classification: {abbrev} (shortest_path_between_frags: {len(shortest_path_between_frags) - 2} > {cutoff})")
                return {'category': category, 'type': 'false_positive', 'subtype': 'cutoff', 'abbrev': abbrev}
            intermediary_rings = self.get_rings_between_two_fragments(mol, aidxf1, aidxf2)
            if len(intermediary_rings) == 0:
                abbrev = 'cmo'
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': 'monopodal', 'subtype': '', 'abbrev': abbrev}
            else:
                sssr = [set(x) for x in mol.GetRingInfo().AtomRings()]  # smallest sets of smallest rings
                intermediary_rings = self._filter_intermediary_rings(mol, intermediary_rings, sssr)
                logging.debug(f"len(intermediary_rings)= {len(intermediary_rings)}")
                if len(intermediary_rings) == 1:
                    type = 'bipodal'  # 1 intermediary ring, which are defined by intersection of aidx, so at least always 1 for each fragment!
                    return self._get_fcc_subtype(category, type, aidxf1, aidxf2, intermediary_rings)
                elif len(intermediary_rings) == 2:
                    type = 'tripodal'
                    return self._get_fcc_subtype(category, type, aidxf1, aidxf2, intermediary_rings)
                else:
                    type = 'other'
                    return self._get_fcc_subtype(category, type, aidxf1, aidxf2, intermediary_rings)

    def _filter_intermediary_rings(self, mol: Mol, intermediary_rings: list, sssr: list):
        """Filter the intermediary rings found within a molecule using the Smallest
        Set of Smallest Rings (SSSR). The idea is that if two fragments have 2 intermediary
        rings that are almost identical but for a few atoms, and these atoms actually are
        contained within the same ring, then only one intermediary ring should be counted.
        This will lower the amount of bipodal that are being identified as tripodal and
        tripodal that are identified as unknown connections.

        :param mol: a molecule
        :param intermediary_rings: the intermediary rings of a molecule
        :param sssr: the smallest set of smallest rings of a molecule
        :return: the filtered intermediary rings of a molecule
        """
        # get a dict of all intermediary rings (ir) with their length as keys
        d = {}
        for i, x in enumerate(intermediary_rings):
            ir_id = f"IR_{str(i).zfill(3)}"
            if len(x) in d.keys():
                # attribute an id to each ir for tracking down
                d[len(x)].append((x, ir_id))
                logging.debug(f"{ir_id}: {x}")
            else:
                d[len(x)] = [(x, ir_id)]
                logging.debug(f"{ir_id}: {x}")

        # filter the dict so we consider only values with more than one ir of the same size
        ir_to_check = {}  # ir with common lengths
        for key, val in d.items():
            if len(val) > 1:
                ir_to_check[key] = d[key]

        # exit this functions if no ir to check
        if not ir_to_check:
            return intermediary_rings

        # continue to investigate IR

        # display IR to checks
        for k in ir_to_check.keys():
            # display info
            logging.debug(f"IR to check with size={k}: {', '.join([x[1] for x in ir_to_check[k]])}")
        # attribute id to each SSSR
        for i in range(len(sssr)):
            sssr[i] = (sssr[i], f"SSSR_{str(i).zfill(3)}")
            logging.debug(f"{sssr[i][1]}: {sssr[i][0]}")

        # check each IR of a ggiven n by checking if the atoms that vary
        to_remove = []
        for k in ir_to_check.keys():
            to_remove_curr = []
            for i in range(len(ir_to_check[k])):
                for j in range(i+1, len(ir_to_check[k])):
                    logging.debug(f"Comparing {ir_to_check[k][i][1]} and {ir_to_check[k][j][1]}")
                    diff = ir_to_check[k][i][0] - ir_to_check[k][j][0]
                    [diff.add(x) for x in ir_to_check[k][j][0] - ir_to_check[k][i][0]]
                    logging.debug(f"Variant atom indices: {diff}")
                    # get a dict with idx: atom so we can look up variant atoms neighbors and find out which are connected to each others
                    atoms = {x: mol.GetAtomWithIdx(x) for x in diff}
                    # get the indices of neighbors for every variant atom
                    neighbors = {idx: [x.GetIdx() for x in a.GetNeighbors()] for idx, a in atoms.items()}
                    # get pairwise combinations of variant atoms according to their neighbors
                    combinations = itertools.combinations([(idx1, idx2) for idx1, idx2 in itertools.combinations(neighbors, 2) if idx2 in neighbors[idx1]], 2)
                    # flatten the subtuple so we can consider combinations for every rings
                    combinations = [tuple(itertools.chain.from_iterable(c)) for c in combinations]
                    [logging.debug(f"Possible combination: {c}") for c in combinations]
                    # check what combinations are found within a ring of the molecule
                    combinations_identified = []
                    for c in combinations:
                        combinations_identified += [c for r in sssr if set(c).issubset(r[0])]
                    if combinations_identified == combinations:
                        logging.debug(f"All combinations were identified, {ir_to_check[k][i][1]} and {ir_to_check[k][j][1]} are equivalent")
                        to_remove_curr.append(ir_to_check[k][i])
                        to_remove_curr.append(ir_to_check[k][j])

            to_remove_curr = [(frozenset(tr[0]), tr[1]) for tr in to_remove_curr]
            to_remove_curr = list(set(to_remove_curr))
            to_remove_ids = list(set([tr[1] for tr in to_remove_curr]))
            logging.debug(f"IR to remove for n={k}: {', '.join([tri for tri in to_remove_ids])}")
            # to_remove_curr = [tr[0] for tr in to_remove_curr]  # get rid of the ids
            # in case all ir of this size are equivalent, just retrieve the first one
            if len(to_remove_curr) == len(ir_to_check[k]):
                logging.debug(f"All IR were detected equivalent, so keeping {to_remove_curr[0][1]}")
                to_remove_curr.pop(0)
                # add current k to the whole mask
                to_remove += to_remove_curr
        to_remove_ids = list(set([tr[1] for tr in to_remove]))
        logging.debug(f"Total IR to remove: {', '.join([tri for tri in to_remove_ids])}")
        # clear to_remove from ids for easier comparison  (maybe lambda funct would perform better here?)
        to_remove = [set(tr[0]) for tr in to_remove]
        # filter the IR by to_remove
        remaining_ir = [ir for ir in intermediary_rings if frozenset(ir) not in to_remove]
        logging.debug(f"Number of remaining_ir: {len(remaining_ir)}")  # we don't have the ids here
        # return the filtered IR
        return remaining_ir

    def _get_fcc_subtype(self, category: str, type: str, aidxf1: set, aidxf2: set, intermediary_rings: list) -> tuple:
        """Return the subtype (spiro, edge, bridged) for bipodal, tripodal and unknown connections.

        :param category: the fragment combination category (connection)
        :param type: the fragment combination type (bipodal, tripodal or unknown)
        :param aidxf1: the atom indices of the first fragment found in the molecule
        :param aidxf2: the atom indices of the second fragment found in the molecule
        :param intermediary_rings: the list of intermediary rings between both fragments and defined by atom indices
        :return: the dictionary specifying fragment combination category, type, subtype and abbrev
        """
        logging.debug("Iterating over remaining IR (new ids and will continue iteration only if subtype=edge):")
        ir_ids = [f"IR_{str(i).zfill(3)}" for i in range(len(intermediary_rings))]
        for i, ir in enumerate(intermediary_rings):
            logging.debug(f"{ir_ids[i]}: {intermediary_rings[i]}")
            intersect_1 = ir.intersection(aidxf1)
            intersect_2 = ir.intersection(aidxf2)
            abbrev = category[0] + type[0]
            logging.debug(f"intersect_1: {intersect_1} ({len(intersect_1)}), intersect_2: {intersect_2} ({len(intersect_2)})")
            if len(intersect_1) == 1 or len(intersect_2) == 1:
                abbrev += 's'
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': type, 'subtype': 'spiro', 'abbrev': abbrev}
            elif 3 <= len(intersect_1) <= 5 or 3 <= len(intersect_2) <= 5:
                abbrev += 'b'
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': type, 'subtype': 'bridged', 'abbrev': abbrev}
            elif len(intersect_1) > 5 or len(intersect_2) > 5:
                abbrev += 'o'  # other
                logging.debug(f"Classification: {abbrev}")
                return {'category': category, 'type': type, 'subtype': 'other', 'abbrev': abbrev}
        # edge otherwise
        abbrev += 'e'
        logging.debug(f"Classification: {abbrev}")
        return {'category': category, 'type': type, 'subtype': 'edge', 'abbrev': abbrev}

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
        df_aidxf['idxf'] = df_aidxf.groupby(['idm', 'idf', 'aidxf_str']).grouper.group_info[0]  # add a seq number that get increased for every group
        # logging.info(f"\n\ndf_aidxf['idxf']:\n {df_aidxf['idxf']}")  # !!! think about what I really want here. For now I just know I don't want this behavior

        # classify fragment combinations
        for gid, g in df_aidxf.groupby('idm'):
            mol = df_mols.loc[gid]['mol']
            hac = mol.GetNumAtoms()  # so we can estimate how well-covered is the molecule by its fragment combinations
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
                    logging.debug("="*80)
                    logging.debug(f"Classifying m={gid}, f1={idf1}:{idxf1}, f2={idf2}:{idxf2}")
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
                    d_fcc['hac'] = hac
                    ds_fcc.append(d_fcc)
        # dataframe with columns in given order
        return DataFrame(ds_fcc, columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', 'aidxf1', 'aidxf2', 'hac'])

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
        num_fcc_ini = len(df_fcc.index)
        df_fcc = df_fcc[df_fcc['abbrev'] != 'cfc']
        logging.debug(f"Removed {len(df_fcc.index)}/{num_fcc_ini} fragment combinations")

        # drop fragments combinations paired with a substructure
        logging.debug(f"Removing substructures from fragment combinations")
        df_substructures = df_fcc[df_fcc['abbrev'] == 'ffs']  # all the substructures in the whole dataframe
        logging.debug(f"Number of substructures found in df_fcc: {len(df_substructures.index)}/{len(df_fcc.index)}")
        logging.debug(f"Substructure combinations:\n\n{df_substructures[['idm', 'fid1', 'fid2', 'abbrev']]}\n")
        logging.debug(f"Determining what fragments should be removed:")
        if len(df_substructures) > 0:
            fid_to_remove = set()
            for gid, g in df_fcc[df_fcc['idm'].isin(df_substructures['idm'])].groupby('idm'):  # iterate only on the groups with at least one substructure
                for rowid, row in g[g['abbrev'] == 'ffs'].iterrows():
                    if len(row['aidxf1']) > len(row['aidxf2']):
                        fid_to_remove.add(row['fid2'])
                    else:
                        fid_to_remove.add(row['fid1'])
                    logging.debug(f"{gid}: {row['fid1']} - {row['fid2']} ==> to_remove={fid_to_remove}")
            # filter the unwanted fragment combinations
            fid_to_remove = list(fid_to_remove)
            logging.debug(f"Number of fragments to remove: {len(fid_to_remove)}")
            nb_fcc_ini = len(df_fcc.index)
            df_fcc = df_fcc[~df_fcc['fid1'].isin(fid_to_remove)]
            df_fcc = df_fcc[~df_fcc['fid2'].isin(fid_to_remove)]
            logging.debug(f"Number of fragment combinations remaining: {len(df_fcc)}/{nb_fcc_ini}")

        if len(df_fcc.index) == 0:
            logging.warning("No fragment remaining for mapping!")
            return None

        return df_fcc

    def map_frags(self, df_fcc: DataFrame, min_frags: int = 2, max_frags: int = 5, max_overlaps: int = 5) -> DataFrame:
        """
        This method process a fragment combinations computed with classify_fragment_combinations
        and return a new DataFrame with a fragment map for each molecule.

        This fragment map is a single line string representation of the fragment connectivity
        within a molecule and follows following syntax:

            >>> f1:0[abbrev1]f2:0-f1:0[abbrev2]f3:0-f2:0[abbrev3]f3:0

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
            logging.debug(f"Number of overlaps found for molecule {gid}: {noverlaps}")
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
                alt_combinations = [list(x) + common_combinations for x in list(itertools.product(*alt_combinations))]

                dfs_fcc_clean = []
                for alt in alt_combinations:
                    df_alt = g[(g['fid1'].isin(alt)) | (g['fid2'].isin(alt))]
                    to_add = True
                    # check if this fcc is already recorded, if so do not record it again
                    for df_fcc_clean in dfs_fcc_clean:
                        if df_alt.equals(df_fcc_clean):
                            to_add = False
                    if to_add:
                        dfs_fcc_clean.append(df_alt)
            else:
                dfs_fcc_clean = [g]
            # fragment map string representation
            for i, df_fcc_clean in enumerate(dfs_fcc_clean):
                # string representation of the fragment combinations of this map
                frag_map_str = '-'.join(list(df_fcc_clean['fid1'].map(str) + "[" + df_fcc_clean['abbrev'] + "]" + df_fcc_clean['fid2'].map(str)))
                # frags: all occurrences of all fragments (f1:0, f1:1, f2:0, etc.)
                frags = list(set(list(df_fcc_clean['fid1'].map(str).values) + list(df_fcc_clean['fid2'].map(str).values)))
                nfrags = len(frags)
                # frags_u: count only different fragments (f1, f2, etc.)
                frags_u = list(set([x.split(":")[0] for x in frags]))
                nfrags_u = len(frags_u)
                # filter results by min/max number of fragments
                if nfrags < min_frags or nfrags > max_frags:
                    logging.debug(f"{gid}: discarding one fragment map because of unsuitable number of fragments ({nfrags})")
                # combine aidxfs from all fragments
                aidxfs = list(df_fcc_clean['aidxf1'].values) + list(df_fcc_clean['aidxf2'].values)
                # organize aidxfs
                aidxfs = dict(zip(frags, aidxfs))  # aidxfs is now a dict with frag: aidfx
                for k in aidxfs.keys():
                    aidxfs[k] = list(aidxfs[k])
                # compute fragment coverage of the molecule
                hac_mol = g.iloc[0]['hac']  # same hac for all entries of the same molecule
                hac_frags = len(list(set([item for sublist in aidxfs.values() for item in sublist])))
                perc_mol_cov_frags = round((hac_frags / hac_mol), 2) * 100

                # avoid issues with pandas and complex data structures by dumping it as string
                aidxfs = json.dumps(aidxfs)
                comb = list(df_fcc_clean['abbrev'].values)
                ncomb = len(comb)
                comb_u = list(set(comb))
                ncomb_u = len(comb_u)
                ds_map.append({'idm': gid, 'fmid': str(i+1).zfill(3), 'nfrags': nfrags, 'nfrags_u': nfrags_u, 'ncomb': ncomb, 'ncomb_u': ncomb_u, 'hac_mol': hac_mol, 'hac_frags': hac_frags, 'perc_mol_cov_frags': perc_mol_cov_frags,  'frags': frags, 'frags_u': frags_u, 'comb': comb, 'comb_u': comb_u, 'aidxfs': aidxfs, 'map_str': frag_map_str})

        # df_map
        return DataFrame(ds_map, columns=['idm', 'fmid', 'nfrags', 'nfrags_u', 'ncomb', 'ncomb_u', 'hac_mol', 'hac_frags', 'perc_mol_cov_frags', 'frags', 'frags_u', 'comb', 'comb_u', 'aidxfs', 'map_str'])
