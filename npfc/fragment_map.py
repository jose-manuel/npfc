"""
Module fragment_map
===================
This modules contains the functions for mapping fragments combinations.
"""

# standard
import logging
import itertools
# data handling
import pandas as pd
from pandas import DataFrame
from collections import OrderedDict
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import AllChem
# graph
import networkx as nx
# docs
from typing import List
# dev
from npfc import draw


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _get_incompatible_fragments_dict(df_overlaps: DataFrame) -> dict:
    """Compute a dictionary indicating what fragment id is incompatible with another:
    d[frag1] = [frag2, frag3, ...]
    d[frag2] = [frag1, frag3, ...]
    d[frag3] = [frag1, frag2, ...]

    :param df_overlaps: DataFrame with only fragment combinations of type ffo (overlaps.)
    :return: a dictionary for easy check of incompatibilities
    """
    d_incompatible = {}
    for i in range(len(df_overlaps.index)):
        row = df_overlaps.iloc[i]
        fid1 = row['fid1']
        fid2 = row['fid2']

        if fid1 not in d_incompatible.keys():
            d_incompatible[fid1] = [fid2]
        else:
            d_incompatible[fid1].append(fid2)

        if fid2 not in d_incompatible.keys():
            d_incompatible[fid2] = [fid1]
        else:
            d_incompatible[fid2].append(fid1)

    return d_incompatible


def _split_overlaps(df_fc: DataFrame, max_overlaps: int) -> DataFrame:
    """Split a fragment combination DataFrame containing overlap entries into different fragment combination DataFrames.

    This function is used within a loop, hence the need for gid parameter (logging).

    :param df_fc: fragment combination DataFrame
    :max_overlaps: maximum number of authorized overlap entries in the DataFrame, if observed number is higher, then an empty DataFrame is returned
    :return: a List of DataFrames
    """
    # Split df_fc into 3 dfs:
    # 1. overlaps: fc of type ffo
    # 2. commons: fc not involving any fragment in ffo
    # 3. variants: fc involving fragments in ffo
    # overlaps
    df_overlaps = df_fc[df_fc['abbrev'] == 'ffo']
    noverlaps = len(df_overlaps.index)
    logging.debug(f"Number of overlaps found: {noverlaps}'")

    if noverlaps == 0:
        return ([df_fc], noverlaps)
    elif noverlaps > max_overlaps:
        logging.debug(f"Too many overlap combinations ({noverlaps}/{max_overlaps}), discarding molecule.")
        return ([], noverlaps)

    # code below is for cases with at least 1 overlap
    d_incompatible = _get_incompatible_fragments_dict(df_overlaps)
    # common
    df_common = df_fc[(~df_fc['fid1'].isin(d_incompatible.keys())) & (~df_fc['fid2'].isin(d_incompatible.keys()))]
    logging.debug(f"Number of common combinations: {len(df_common.index)}")
    # variants
    df_variants = df_fc[ (df_fc['abbrev'] != 'ffo') & ((df_fc['fid1'].isin(d_incompatible.keys()) | (df_fc['fid2'].isin(d_incompatible.keys()))))]
    logging.debug(f"Number of variant combinations: {len(df_variants.index)}")

    # define a list of sets of alternative fragments (i.e. [(A, B), (C, D)]) for computing all possible alternative fmap possibilities
    alt_frags = list(set([tuple(sorted([k] + v)) for k, v in d_incompatible.items()]))

    # compute 1 df for each alternative fmap and then concatenate it all into one single df
    dfs_alt_curr = []
    for product in itertools.product(*alt_frags):
        logging.debug(f"Current alternative route for fmaps: {product}")
        to_remove = []
        for p in product:
            to_remove += d_incompatible[p]
        df_variants_curr = df_variants[((~df_variants['fid1'].isin(to_remove)) & (~df_variants['fid2'].isin(to_remove)))]
        dfs_alt_curr.append(pd.concat([df_common, df_variants_curr]))

    return (dfs_alt_curr, noverlaps)


def _split_unconnected(dfs_fcc_clean: List[DataFrame]) -> List[DataFrame]:
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


def generate(df_fcc: DataFrame, min_frags: int = 2, max_frags: int = 5, max_overlaps: int = 5) -> DataFrame:
    """This method process a fragment combinations DataFrame
    and return a new DataFrame with a fragment map for each molecule.

    |pic1| |pic2|

    .. |pic1| image:: _images/map_ex1_mol.svg
       :width: 43%

    .. |pic2| image:: _images/map_ex1_graph.svg
       :width: 56%

    Each highlighted fragment of the molecule consists of a node of the graph.
    Fragment Combinations are used as edges for displaying fragment interactions and are annotated
    with the idm as well as the category of the combination.

    A str representation is also computed:
        frag1[cmo]frag2-frag2[fed]frag3

    Two objects are computed and stored as b64 strings:
        - colormap: a custom object with 3 informations regarding highlight colors (RGB values):
            - fragments: the color attributed to each fragment
            - atoms: the color attributed to each atom
            - bonds: the color attributed to each bond
        - graph: a nx object used comparing fragment connectivity among molecules.

    Molecules can be filtered using thresholds.

    :param df_fcc: a Dataframe with pairwise fragment combinations
    :param min_frags: a threshold for the minimum number of fragments allowed per fragment map
    :param max_frags: a threshold for the maximum number of fragments allowed per fragment map
    :param max_overlaps: a threshold for the maximum number of overlap combinations found in the molecule
    :return: a DataFrame representing fragment maps
    """
    # split by overlaps

    logging.debug(f"Mapping fragments")

    ds_map = []
    for gid, g in df_fcc.groupby('idm'):
        logging.debug(f"Current Molecule: {gid}")
        # split overlaps into different Dataframes
        dfs_fcc_clean, noverlaps = _split_overlaps(g, max_overlaps)
        logging.debug(f"dfs_fcc_clean\n{dfs_fcc_clean}\n")
        if len(dfs_fcc_clean) == 0:
            continue

        # compute fragment connectivity graph objects so we can split up disconnected subgraphs
        dfs_fcc_ready = _split_unconnected(dfs_fcc_clean)
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
                    d_aidxs[row["idf1"]] = [row["_aidxf1"]]
                elif row["_aidxf1"] not in d_aidxs[row["idf1"]]:
                    d_aidxs[row["idf1"]].append(row["_aidxf1"])
                # idf2
                if row["idf2"] not in d_aidxs.keys():
                    d_aidxs[row["idf2"]] = [row["_aidxf2"]]
                elif row["_aidxf2"] not in d_aidxs[row["idf2"]]:
                    d_aidxs[row["idf2"]].append(row["_aidxf2"])

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
            edge_attr = ['abbrev', 'n_abbrev', 'idm', 'n_active', 'n_tot', 'success_rate']
            edge_attr = [x for x in edge_attr if x in df_fcc_clean.columns]
            graph = nx.from_pandas_edgelist(df_fcc_clean, source="idf1", target="idf2", edge_attr=edge_attr)

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
            ds_map.append({'idm': gid, 'fmid': str(i+1).zfill(3), 'nfrags': nfrags, 'nfrags_u': nfrags_u, 'ncomb': ncomb, 'ncomb_u': ncomb_u, 'hac_mol': hac_mol, 'hac_frags': hac_frags, 'perc_mol_cov_frags': perc_mol_cov_frags, 'frags': frags, 'frags_u': frags_u, 'comb': comb, 'comb_u': comb_u, 'fmap_str': frag_map_str, '_d_aidxs': d_aidxs, '_colormap': colormap, '_fmap': graph, 'mol': mol})

    # df_map
    return DataFrame(ds_map, columns=['idm', 'fmid', 'nfrags', 'nfrags_u', 'ncomb', 'ncomb_u', 'hac_mol', 'hac_frags', 'perc_mol_cov_frags', 'frags', 'frags_u', 'comb', 'comb_u', 'fmap_str', '_d_aidxs', '_colormap', '_fmap', 'mol'])