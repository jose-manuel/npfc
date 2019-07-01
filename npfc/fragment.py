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
# data handling
from collections import OrderedDict
from collections import Counter
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import AllChem
# graph
import networkx as nx
# docs
from pandas import DataFrame
from typing import List
# dev
from npfc import draw

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Matcher:
    """Create a Matcher object."""

    def __init__(self):
        pass

    def run(self,
            df_mols: DataFrame,
            df_frags: DataFrame,
            col_mol_mols: str = 'mol',
            col_mol_frags: str = 'mol',
            ) -> DataFrame:
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
        d['idxf'] = []
        d['aidxf'] = []
        d['mol_perc'] = []  # proportion of the molecule the substructure represents
        d['mol'] = []  # encode the molecule here so we don't have to combine multiple files when trying to have a look at the results
        d['mol_frag'] = []  # strucutre of the fragment
        # begin
        for idm, rowm in df_mols.iterrows():
            mol = Mol(rowm[col_mol_mols])
            hac = rowm[col_mol_mols].GetNumAtoms()
            for idf, rowq in df_frags.iterrows():
                # perform the substructure search on mol so the latest matching fragment does not get highlighted
                matches = mol.GetSubstructMatches(rowq[col_mol_frags])
                for i, m in enumerate(matches):
                    d['idm'].append(idm)
                    d['idf'].append(idf)
                    d['aidxf'].append(frozenset(m))  # frozenset so we can use intersection, etc. and still remove dupl. easily
                    d['idxf'].append(str(i))
                    d['mol_perc'].append(round(len(m)/hac, 2) * 100)
                    d['mol'].append(rowm[col_mol_mols])
                    d['mol_frag'].append(rowq[col_mol_frags])

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
                                       df_aidxf: DataFrame,
                                       cutoff: int = 3,
                                       clean: bool = True) -> DataFrame:
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

        :param df_aidxf: the input DataFrame with substructure matches
        :param cutoff: the maximum number of intermediary atoms between 2 fragments
        :param clean: remove false positives such as cfc or substructures from results by calling the clean method.
        :return: a DataFrame with all fragment combination classifications
        """
        ds_fcc = []
        logging.debug(df_aidxf.columns)
        # labelling idxf
        df_aidxf['aidxf_str'] = df_aidxf['aidxf'].map(str)  # sets are an unhashable type...

        # logging.info(f"\n\ndf_aidxf['idxf']:\n {df_aidxf['idxf']}")  # !!! think about what I really want here. For now I just know I don't want this behavior

        # classify fragment combinations
        for gid, g in df_aidxf.groupby('idm'):
            mol = g.iloc[0]['mol']
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
                    d_fcc['mol'] = mol
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
        logging.debug("="*80)
        # dataframe with columns in given order
        df_fcc = DataFrame(ds_fcc, columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', 'aidxf1', 'aidxf2', 'hac', 'mol'])
        # clean results from false positives
        if clean:
            return self.clean(df_fcc)
        return df_fcc

    def clean(self, df_fcc: DataFrame) -> DataFrame:
        """Clean a df_fcc by removing false positives such as substructures and
        cutoff combinations.

        Overlaps are still kept for further processing.

        :param df_fcc: a fcc DataFrame
        :return: a cleaned fcc DataFrame
        """

        # clean the data

        logging.debug("Now cleaning fragment combinations")

        # drop cutoff combinations
        logging.debug(f"Removing cutoff connections from fragment combinations")
        num_fcc_ini = len(df_fcc.index)
        df_fcc = df_fcc[df_fcc['abbrev'] != 'cfc']
        logging.debug(f"Number of remaining fragment combinations: {len(df_fcc.index)}/{num_fcc_ini}")

        # drop fragments combinations paired with a substructure
        logging.debug(f"Removing substructures from fragment combinations")
        df_substructures = df_fcc[df_fcc['abbrev'] == 'ffs']  # all the substructures in the whole dataframe
        num_substructures = len(df_substructures.index)
        logging.debug(f"Number of substructures found in df_fcc: {num_substructures}/{len(df_fcc.index)}")
        # in case of substructures to remove, iterate over all identified subtructures for each molecule,
        # determine what fragments are part of others and discard all entries with them
        if num_substructures > 0:
            logging.debug(f"Substructure combinations:\n\n{df_substructures[['idm', 'fid1', 'fid2', 'abbrev']]}\n")
            logging.debug(f"Determining what fragments should be removed:")
            # intialize the iteration
            rowids_to_remove = []  # the rowids of the df_fcc dataframe to remove
            for gid, g in df_fcc[df_fcc['idm'].isin(df_substructures['idm'])].groupby('idm'):  # iterate only on the groups with at least one substructure
                fid_to_remove = set()   # fid of substructures identified for the current molecule
                # for each molecule, look at what fids we should remove
                for rowid, row in g[g['abbrev'] == 'ffs'].iterrows():
                    # combination ifs ffs, so remove either fid1 or fid2 depending on hac
                    if len(row['aidxf1']) > len(row['aidxf2']):
                        fid_to_remove.add(row['fid2'])
                    else:
                        fid_to_remove.add(row['fid1'])
                    # display some debugging
                    logging.debug(f"{gid}: {row['fid1']} - {row['fid2']} ==> fid_to_remove={fid_to_remove}")
                    # register df_fcc rowids that will be removed for this substructure
                    rowids_to_remove += list(g[g["fid1"].isin(list(fid_to_remove))].index) + list(g[g["fid2"].isin(list(fid_to_remove))].index)
            # remove dupl in rowids_to_remove
            rowids_to_remove = list(set(rowids_to_remove))
            # filter the unwanted fragment combinations
            logging.debug(f"Number of fragments combinations to remove: {len(rowids_to_remove)}")
            nb_fcc_ini = len(df_fcc.index)
            df_fcc = df_fcc.loc[~df_fcc.index.isin(rowids_to_remove)]
            logging.debug(f"Number of fragment combinations remaining: {len(df_fcc)}/{nb_fcc_ini}")

        return df_fcc

    def _split_overlaps(self, gid: str, g: DataFrame, max_overlaps: int) -> DataFrame:
        """
        Split a fcc DataFrame containing overlap entries into different DataFrames.

        This function is used within a loop, hence the need for gid parameter (logging).

        :param gid: current group id
        :param g: current group DataFrame
        :max_overlaps: maximum number of authorized overlap entries in the DataFrame, if observed number is higher, then an empty DataFrame is returned
        :return: a List of DataFrames
        """
        # entries with an overlap
        overlaps = g[g['abbrev'] == 'ffo']
        noverlaps = len(overlaps.index)
        if noverlaps > 0:  # the code below could certainly be improved, but this case should not happen too often
            logging.debug(f"Number of overlaps found for molecule {gid}: {noverlaps}")
            # filter out molecules with too many overlaps
            if noverlaps > max_overlaps:
                logging.debug(f"Too many overlap combinations ({noverlaps}), discarding molecule '{gid}'")
                # return empty DataFrameas well as the number of overlaps found
                return (DataFrame(columns=g.columns), noverlaps)
            # remove overlaps from the current group
            g = g[g['abbrev'] != 'ffo']
            # get fragment ids of the common parts of alternative paths
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
                # check if this fcc is already recorded, if so do not record it again,
                # this is useful in case of common parts being the only remaining parts of overlaps
                for df_fcc_clean in dfs_fcc_clean:
                    if df_alt.equals(df_fcc_clean):
                        to_add = False
                if to_add:
                    dfs_fcc_clean.append(df_alt)
        else:
            dfs_fcc_clean = [g]

        return (dfs_fcc_clean, noverlaps)

    def _split_unconnected(self, dfs_fcc_clean: List[DataFrame]) -> List[DataFrame]:
        """
        From a list of DataFrames of Fragment Combinations, split up DataFrames containing
        two or more unconnected parts in different DataFrames.
        For instance if a molecule has two paired fragments left and two fragments right
        too far from other, then we get something like 2 and 2 instead of 4.

        :param dfs_fcc_clean: a List of fcc DataFrames
        :return: an updated List of connected fcc DataFrames
        """
        dfs_fcc_ready = []
        for i, df_fcc_clean in enumerate(dfs_fcc_clean):
            # compute a graph with each fid as a node, one row means an edge between 2 fid
            fc_graph = nx.from_pandas_edgelist(df_fcc_clean, "fid1", "fid2")
            fc_subgraphs = list(nx.connected_component_subgraphs(fc_graph))
            num_fc_subgraphs = len(fc_subgraphs)
            # splitting up subgraphs
            if num_fc_subgraphs > 1:
                logging.debug(f"Fragment Connectivity" + f"{i}".rjust(5) + f": found {num_fc_subgraphs} fc_subgraphs, so splitting up")
                # for each subgraph, record corresponding rows in df only
                for fc_subgraph in fc_subgraphs:
                    nodes = list(fc_subgraph.nodes())
                    df_fcc_subgraph = df_fcc_clean[((df_fcc_clean['fid1'].isin(nodes)) | (df_fcc_clean['fid2'].isin(nodes)))].copy()
                    dfs_fcc_ready.append(df_fcc_subgraph)
            else:
                dfs_fcc_ready.append(df_fcc_clean)

        return dfs_fcc_ready

    def map_frags(self, df_fcc: DataFrame, min_frags: int = 2, max_frags: int = 5, max_overlaps: int = 5) -> DataFrame:
        """This method process a fragment combinations computed with classify_fragment_combinations
        and return a new DataFrame with a fragment map for each molecule.

        This fragment map is a single line string representation of the fragment connectivity
        within a molecule and follows following syntax:

            >>> f1:0[cmo]f2:0-f1:0[fed]f3:0-f2:0[cmo]f3:0

        with a fragment id being composed of 2 parts seperated by ":":
            - f1: the fragment type
            - 0: the occurrence number of the fragment type in this molecule

        """
        # split by overlaps

        logging.debug(f"Mapping fragments")

        ds_map = []
        for gid, g in df_fcc.groupby('idm'):

            # split overlaps into different Dataframes
            dfs_fcc_clean, noverlaps = self._split_overlaps(gid, g, max_overlaps)
            if len(dfs_fcc_clean) == 0:
                continue

            # compute fragment connectivity graph objects so we can split up disconnected subgraphs
            dfs_fcc_ready = self._split_unconnected(dfs_fcc_clean)

            # compute the entries of the df_map
            for i, df_fcc_clean in enumerate(dfs_fcc_ready):

                # string representation of the fragment combinations of this map
                frag_map_str = '-'.join(list(df_fcc_clean['fid1'].map(str) + "[" + df_fcc_clean['abbrev'] + "]" + df_fcc_clean['fid2'].map(str)))

                # d_aidxs: a dict containing the occurrences of each fragment type
                d_aidxs = {}  # normal dict
                for j in range(len(df_fcc_clean.index)):
                    row = df_fcc_clean.iloc[j]
                    # idf1
                    if row["idf1"] not in d_aidxs.keys():
                        d_aidxs[row["idf1"]] = [row["aidxf1"]]
                    elif row["aidxf1"] not in d_aidxs[row["idf1"]]:
                        d_aidxs[row["idf1"]].append(row["aidxf1"])
                    # idf2
                    if row["idf2"] not in d_aidxs.keys():
                        d_aidxs[row["idf2"]] = [row["aidxf2"]]
                    elif row["aidxf2"] not in d_aidxs[row["idf2"]]:
                        d_aidxs[row["idf2"]].append(row["aidxf2"])

                # sort d_aidxs for reproducible colormaps
                d_aidxs = OrderedDict(sorted(d_aidxs.items()))

                # count fragment occurrences (non-unique)
                frags = list(set([x for x in df_fcc_clean['fid1'].map(str).values] + [x for x in df_fcc_clean['fid2'].map(str).values]))
                nfrags = len(frags)

                # filter results by min/max number of fragment occurrences
                if nfrags < min_frags or nfrags > max_frags:
                    logging.debug(f"{gid}: discarding one fragment map because of unsuitable number of fragments ({nfrags})")
                    continue

                # count unique fragment types (unique)
                frags_u = list(d_aidxs.keys())
                nfrags_u = len(frags_u)

                # compute fragment coverage of the molecule
                hac_mol = g.iloc[0]['hac']  # same hac for all entries since this is the same molecule anyway
                # hac
                hac_frags = set()
                for k in d_aidxs.keys():
                    hac_frags.update(set([item for sublist in d_aidxs[k] for item in sublist]))
                hac_frags = len(hac_frags)
                # perc
                perc_mol_cov_frags = round((hac_frags / hac_mol), 2) * 100

                # compute a new graph again but this time on a single subgraph and with edge labels (room for optimization)
                # count the number of equivalent edges (sames ids and same abbrev)
                df_fcc_clean = df_fcc_clean.copy()  # ### one day I will have to understand why all of the Pandas warnings appear all the time
                df_fcc_clean['n_abbrev'] = df_fcc_clean.groupby(['idf1', 'idf2', 'abbrev'])['abbrev'].transform('count')
                df_fcc_clean.drop_duplicates(subset=["idf1", "idf2", "abbrev"], keep="first", inplace=True)
                # compute the graph
                graph = nx.from_pandas_edgelist(df_fcc_clean, source="idf1", target="idf2", edge_attr=["abbrev", "n_abbrev"])

                # same molecule in each row, so to use the first one is perfectly fine
                mol = df_fcc_clean.iloc[0]['mol']

                # in case of overlaps, the same molecule will be used more than once,
                # so make a copy of the original so highlights are truly independant
                if noverlaps > 0:
                    mol = Mol(mol)

                # attribute colors to each fragment atoms/bonds
                colormap = draw.ColorMap(mol, d_aidxs, draw.colors)

                comb = list(df_fcc_clean['abbrev'].values)
                ncomb = len(comb)
                comb_u = list(set(comb))
                ncomb_u = len(comb_u)
                ds_map.append({'idm': gid, 'fmid': str(i+1).zfill(3), 'nfrags': nfrags, 'nfrags_u': nfrags_u, 'ncomb': ncomb, 'ncomb_u': ncomb_u, 'hac_mol': hac_mol, 'hac_frags': hac_frags, 'perc_mol_cov_frags': perc_mol_cov_frags, 'frags': frags, 'frags_u': frags_u, 'comb': comb, 'comb_u': comb_u, 'map_str': frag_map_str, 'd_aidxs': d_aidxs, 'colormap': colormap, 'graph': graph, 'mol': mol})

        # df_map
        return DataFrame(ds_map, columns=['idm', 'fmid', 'nfrags', 'nfrags_u', 'ncomb', 'ncomb_u', 'hac_mol', 'hac_frags', 'perc_mol_cov_frags', 'frags', 'frags_u', 'comb', 'comb_u', 'map_str', 'd_aidxs', 'colormap', 'graph', 'mol'])
