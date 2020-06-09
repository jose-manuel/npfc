"""
Module fragment_combination_graph
===================================
This modules contains the functions for generating fragments combination graphs
for individual molecules.
"""

# standard
import logging
import itertools
import toolz
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
from typing import Union
# dev
from npfc import draw
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# for computing fragment graphs
DF_FG_COLS = ['idm',
              'idfcg',
              'nfrags',
              'nfrags_u',
              'ncomb',
              'ncomb_u',
              'hac_mol',
              'hac_frags',
              'perc_mol_cov_frags',
              '_frags',
              '_frags_u',
              'comb',
              'comb_u',
              'fcg_str',
              '_d_aidxs',
              '_colormap',
              '_fcg',
              'mol',
              '_d_mol_frags',
              ]

# for annotating fragment graphs with PNP
DF_PNP_COLS = DF_FG_COLS + ['pnp_fm', 'pnp_mol']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _clear_ffs(df_fcc: DataFrame) -> DataFrame:
    """Clear ffs combinations by discarding any combination in which  the smaller
    fragments is involved.

    :param df_fcc: a fcc DataFrame
    :return: a cleaned fcc DataFrame
    """

    # clean the data

    logging.debug("Now cleaning fragment combinations")

    # drop cutoff combinations
    logging.debug("Removing cutoff connections from fragment combinations")
    num_fcc_ini = len(df_fcc.index)
    logging.debug("Number of remaining fragment combinations: %s/%s", len(df_fcc.index), num_fcc_ini)

    # drop fragments combinations paired with a substructure
    logging.debug("Removing substructures from fragment combinations")
    df_substructures = df_fcc[df_fcc['fcc'] == 'ffs']  # all the substructures in the whole dataframe
    num_substructures = len(df_substructures.index)
    logging.debug("Number of substructures found in df_fcc: %s/%s", num_substructures, len(df_fcc.index))
    # in case of substructures to remove, iterate over all identified subtructures for each molecule,
    # determine what fragments are part of others and discard all entries with them
    if num_substructures > 0:
        logging.debug("Substructure combinations:\n\n%s\n", df_substructures[['idm', 'fc']])
        logging.debug("Determining what fragments should be removed:")
        # intialize the iteration
        rowids_to_remove = []  # the rowids of the df_fcc dataframe to remove
        for gid, g in df_fcc[df_fcc['idm'].isin(df_substructures['idm'])].groupby('idm'):  # iterate only on the groups with at least one substructure
            fid_to_remove = set()   # fid of substructures identified for the current molecule
            # for each molecule, look at what fids we should remove
            for rowid, row in g[g['fcc'] == 'ffs'].iterrows():
                # combination ifs ffs, so remove either fid1 or fid2 depending on hac
                if len(row['_aidxf1']) > len(row['_aidxf2']):
                    fid_to_remove.add(row['fid2'])
                else:
                    fid_to_remove.add(row['fid1'])
                # display some debugging
                logging.debug("%s: %s - %s ==> fid_to_remove: %s", gid, row['fid1'], row['fid2'], fid_to_remove)
                # register df_fcc rowids that will be removed for this substructure
                rowids_to_remove += list(g[g["fid1"].isin(list(fid_to_remove))].index) + list(g[g["fid2"].isin(list(fid_to_remove))].index)
        # remove dupl in rowids_to_remove
        rowids_to_remove = list(set(rowids_to_remove))
        # filter the unwanted fragment combinations
        logging.debug("Number of fragments combinations to remove: %s", len(rowids_to_remove))
        nb_fcc_ini = len(df_fcc.index)
        df_fcc = df_fcc.loc[~df_fcc.index.isin(rowids_to_remove)]
        logging.debug("Number of fragment combinations remaining: %s/%s", len(df_fcc), nb_fcc_ini)

    return df_fcc


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
    df_overlaps = df_fc[df_fc['fcc'] == 'ffo']
    noverlaps = len(df_overlaps.index)
    logging.debug("Number of overlaps found: %s'", noverlaps)

    if noverlaps == 0:
        return ([df_fc], noverlaps)
    elif noverlaps > max_overlaps:
        logging.debug("Too many overlap combinations (%s/%s), discarding molecule.", noverlaps, max_overlaps)
        return ([], noverlaps)

    # code below is for cases with at least 1 overlap
    d_incompatible = _get_incompatible_fragments_dict(df_overlaps)
    # common
    df_common = df_fc[(~df_fc['fid1'].isin(d_incompatible.keys())) & (~df_fc['fid2'].isin(d_incompatible.keys()))]
    logging.debug("Identified %s variant combinations:\n%s", len(df_common.index), df_common[['idm', 'fc']])

    # variants
    df_variants = df_fc[(df_fc['fcc'] != 'ffo') & ((df_fc['fid1'].isin(d_incompatible.keys()) | (df_fc['fid2'].isin(d_incompatible.keys()))))]
    logging.debug("Identified %s variant combinations:\n%s", len(df_variants.index), df_variants[['idm', 'fcc']])

    # define a list of sets of alternative fragments (i.e. [(A, B), (C, D)]) for computing all possible alternative fgraph possibilities
    alt_frags = list(set([tuple(sorted([k] + v)) for k, v in d_incompatible.items()]))
    # remove sublists included in others: i.e. ('1724:0', '627:1'), ('1724:0', '627:0'), ('1724:0', '627:0', '627:1') => ('1724:0', '627:0', '627:1')
    alt_frags = [set(x) for x in alt_frags]  # now the former set of tuples is a list of sets
    alt_frags.sort(key=len)  # do only 1 check: left in right and not right in left as well
    alt_frags = [tuple(sorted(list(x))) for x in list(filter(lambda f: not any(f < g for g in alt_frags), alt_frags))]  # discard any set that is subset of another
    logging.debug("Alternative fragments:\n" + '\n'.join([str(x) for x in alt_frags]))

    # compute 1 df for each alternative fgraph and then concatenate it all into one single df
    dfs_alt_curr = []
    for product in itertools.product(*alt_frags):
        logging.debug("Current alternative route for fgraphs: %s", product)
        to_remove = []
        for p in product:
            to_remove += d_incompatible[p]
        df_variants_curr = df_variants[((~df_variants['fid1'].isin(to_remove)) & (~df_variants['fid2'].isin(to_remove)))]
        dfs_alt_curr.append(pd.concat([df_common, df_variants_curr]))

    # unfortunately, there are a lot of duplicates going on after splitting by overlaps
    dfs_alt_curr = list(toolz.unique(dfs_alt_curr, key=lambda x: utils.encode_object(x)))

    # clear duplicate dfs
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
        G = nx.from_pandas_edgelist(df_fcc_clean, "fid1", "fid2")
        fc_subgraphs = list(G.subgraph(c) for c in nx.connected_components(G))
        num_fc_subgraphs = len(fc_subgraphs)
        # splitting up subgraphs
        if num_fc_subgraphs > 1:
            logging.debug("Fragment Connectivity -- %s: found %s subgraphs, so splitting up", i, num_fc_subgraphs)
            # for each subgraph, record corresponding rows in df only
            for fc_subgraph in fc_subgraphs:
                nodes = list(fc_subgraph.nodes())
                df_fcc_subgraph = df_fcc_clean[((df_fcc_clean['fid1'].isin(nodes)) | (df_fcc_clean['fid2'].isin(nodes)))].copy()
                dfs_fcc_ready.append(df_fcc_subgraph)
        else:
            dfs_fcc_ready.append(df_fcc_clean)

    return dfs_fcc_ready


def generate(df_fcc: DataFrame, min_frags: int = 2, max_frags: int = 5, max_overlaps: int = 5, split_unconnected: bool = True, clear_ffs: bool = True, palette: str = None) -> DataFrame:
    """This method process a fragment combinations DataFrame
    and return a new DataFrame with a fragment combination graph for each molecule.

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
    :param min_frags: a threshold for the minimum number of fragments allowed per fragment combination graph
    :param max_frags: a threshold for the maximum number of fragments allowed per fragment combination graph
    :param max_overlaps: a threshold for the maximum number of overlap combinations found in the molecule
    :return: a DataFrame representing fragment combination graphs
    """
    # split by overlaps

    logging.debug("Mapping fragments")

    ds_fcg = []
    for gid, g in df_fcc.groupby('idm'):
        logging.debug("Current Molecule: %s", gid)
        # split overlaps into different Dataframes
        dfs_fcc_clean, noverlaps = _split_overlaps(g, max_overlaps)

        if logging.getLogger().level == logging.DEBUG:
            for i, df_fcc_clean in enumerate(dfs_fcc_clean):
                logging.info("df_fcc_clean #%s\n\n%s\n", i, df_fcc_clean[['idm', 'fc']])

        if len(dfs_fcc_clean) == 0:
            continue

        # compute fragment connectivity graph objects so we can split up disconnected subgraphs
        if split_unconnected:
            dfs_fcc_ready = _split_unconnected(dfs_fcc_clean)
        else:
            dfs_fcc_ready = dfs_fcc_clean

        # clear ffs
        if clear_ffs:
            dfs_fcc_ready = [_clear_ffs(df_fcc_ready) for df_fcc_ready in dfs_fcc_ready]

        # compute the entries of the df_fcg
        for i, df_fcc_clean in enumerate(dfs_fcc_ready):

            # string representation of the fragment combinations of this map
            fragment_combination_graph_str = '-'.join(list(df_fcc_clean['fc'].map(str)))
            logging.debug('fcg_str: %s', fragment_combination_graph_str)

            # d_aidxs: a dict containing the occurrences of each fragment type
            d_aidxs = {}
            d_frags = {}
            for j in range(len(df_fcc_clean.index)):
                row = df_fcc_clean.iloc[j]
                # idf1
                if row["idf1"] not in d_aidxs.keys():
                    d_aidxs[row["idf1"]] = [row["_aidxf1"]]
                    d_frags[row["idf1"]] = row["mol_frag_1"]
                elif row["_aidxf1"] not in d_aidxs[row["idf1"]]:
                    d_aidxs[row["idf1"]].append(row["_aidxf1"])
                # idf2
                if row["idf2"] not in d_aidxs.keys():
                    d_aidxs[row["idf2"]] = [row["_aidxf2"]]
                    d_frags[row["idf2"]] = row["mol_frag_2"]
                elif row["_aidxf2"] not in d_aidxs[row["idf2"]]:
                    d_aidxs[row["idf2"]].append(row["_aidxf2"])

            # sort d_aidxs for reproducible colormaps  ### might be the cause of the wrong coloration in alternative fgraphs due to overlaps
            d_aidxs = OrderedDict(sorted(d_aidxs.items()))
            # count fragment occurrences (non-unique)
            frags = list(set([x for x in df_fcc_clean['fid1'].map(str).values] + [x for x in df_fcc_clean['fid2'].map(str).values]))
            nfrags = len(frags)

            # filter results by min/max number of fragment occurrences
            if nfrags < min_frags or nfrags > max_frags:
                logging.debug("%s: discarding one fragment combination graph because of unsuitable number of fragments (%s)", gid, nfrags)
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
            # count the number of equivalent edges (sames ids and same fcc)
            df_fcc_clean = df_fcc_clean.copy()  # ### one day I will have to understand why all of the Pandas warnings appear all the time
            # here I used to compress all combinations in common
            # df_fcc_clean['n_fcc'] = df_fcc_clean.groupby(['idf1', 'idf2', 'fcc'])['fcc'].transform('count')
            # df_fcc_clean.drop_duplicates(subset=["idf1", "idf2", "fcc"], keep="first", inplace=True)
            # create two new columns: cps and cpt for indicating connection points of source and target
            df_fcc_clean['cps'], df_fcc_clean['cpt'] = zip(*df_fcc_clean['fc'].map(lambda x: (x.split('[')[0].split('@')[1], x.split(']')[1].split('@')[1])))
            # compute the graph
            edge_attr = ['fcc', 'n_fcc', 'idm', 'cps', 'cpt']
            edge_attr = [x for x in edge_attr if x in df_fcc_clean.columns]
            df_fcc_clean = df_fcc_clean.sort_values(['idf1', 'idf2'])
            G = nx.from_pandas_edgelist(df_fcc_clean, source="idf1", target="idf2", edge_attr=edge_attr)
            # same molecule in each row, so to use the first one is perfectly fine
            mol = df_fcc_clean.iloc[0]['mol']

            # in case of overlaps, the same molecule will be used more than once,
            # so make a copy of the original so highlights are truly independant
            if noverlaps > 0:
                mol = Mol(mol)

            # attribute colors to each fragment atoms/bonds
            colormap = draw.ColorMap(mol, d_aidxs, palette)

            comb = list(df_fcc_clean['fcc'].values)
            ncomb = len(comb)
            comb_u = list(set(comb))
            ncomb_u = len(comb_u)
            ds_fcg.append({'idm': gid, 'idfcg': str(i+1).zfill(3), 'nfrags': nfrags, 'nfrags_u': nfrags_u, 'ncomb': ncomb, 'ncomb_u': ncomb_u, 'hac_mol': hac_mol, 'hac_frags': hac_frags, 'perc_mol_cov_frags': perc_mol_cov_frags, '_frags': frags, '_frags_u': frags_u, '_comb': comb, '_comb_u': comb_u, 'fcg_str': fragment_combination_graph_str, '_d_aidxs': d_aidxs, '_colormap': colormap, '_fcg': G, 'mol': mol, '_d_mol_frags': d_frags})

    # put it all together
    df_fcg = DataFrame(ds_fcg, columns=DF_FG_COLS).drop_duplicates(subset=['fcg_str'])
    df_fcg['idfcg'] = df_fcg.groupby('idm').cumcount().map(lambda x: str(x+1).zfill(3))
    # incorporate the idcfg to the graphs
    df_fcg.apply(lambda x: nx.classes.function.set_edge_attributes(x['_fcg'], x['idfcg'], 'idcfg'), axis=1)

    return df_fcg


def has_only_referenced_edges(edges: tuple, edges_ref: tuple) -> bool:
    """Check if at least one edge in edges is not present within edges_ref.
    Edges are tuple of syntax (u, v, d) with d being the dict with the attributes.

    :param edges: the edges of the target molecule fcg
    :param edges_ref: the edges of reference molecule fcg
    :return: False if at least 1 edge is not found in the reference, True otherwise
    """
    for e in edges:
        if e not in edges_ref:
            return False
    return True


def filter_edges_attributes(edges: list, cols: list) -> list:
    """Networkx can either return one or all properties. Now that I am aware of this,
    I extract all attributes and then filter out the ones I do not want to use.

    :param edges: the edges as a list of tuple of syntax (u, v, d) with d being the dict with the attributes
    :param cols: the list of attributes to use for PNP labelling.
    """
    return [(row[0], row[1], {k: v for k, v in row[2].items() if k in cols}) for row in edges]


def get_pnp_references(edges: tuple, df_ref: DataFrame, target_node: frozenset = None) -> tuple:
    """Return a tuple of reference idms. A reference is recorded if all edges of the target molecule are included  at once
    within the reference row edges. In case no reference is found, an empty tuple is returned.

    :param edges: the edges as a list of tuple of syntax (u, v, d) with d being the dict with the attributes
    :param df_ref: the dataframe containing the edges to use for references.
    :param target_nodes: a frozenset of fragment ids found in the target edges to use for filering references to compare (optimization)
    """
    frags_u = frozenset([x[0] for x in edges] + [x[1] for x in edges])

    # avoid unnecessary graph comparisons: no way the syntehtic compound will be matching a NP fc graph if it does not even have all the nodes
    df_ref = df_ref[df_ref['_frags_u'].map(lambda x: frags_u.issubset(x))]

    # get the subset of references that match with this particular molecule
    df_ref = df_ref[df_ref['edges'].map(lambda x: has_only_referenced_edges(edges, x))]

    # if no ref is found, just label as PNP
    if len(df_ref) < 1:
        return tuple()
    else:
        df_ref['idm_idfcg'] = df_ref['idm'].astype(str) + ':' + df_ref['idfcg'].astype(str)
        return tuple(df_ref['idm_idfcg'].values)


def annotate_pnp_fcg(df_fcg, df_fcg_ref, data=['fcc']) -> DataFrame:
    """Search and Identify for PNP molecules in the input DataFrame (df_fcg).
    PNP molecules are defined as molecules containing natural fragments combinations
    that are not found in a reference natural dataset.
    Fragment combinations are defined by extracting information from networkx graphs containing 3 informations:

        - source: id of fragment 1
        - target: id of fragment 2
        - attributes to consider: fcc (by default), cps, cpt, etc.

    Three new columns are appended to the input DataFrame:

        - _pnp_ref: the list of a references found for the target fcg
        - pnp_fcg: True if the input fcg has no match with any reference fcg, False otherwise
        - pnp_mol: True if the input molecule has no matching fcg with any reference fcg, False otherwise

    :param df_fcg: the input DataFrame
    :param df_fcg_ref: the reference DataFrame
    :param data: the list of edge attributes to consider during fcg comparison
    :return: the input DataFrame with 3 new pnp columns
    """

    # extract edge informations for targets
    df_fcg["edges"] = df_fcg["_fcg"].map(lambda x: x.edges(data=True))
    df_fcg["edges"] = df_fcg["edges"].map(lambda x: filter_edges_attributes(x, data))
    df_fcg['_frags_u'] = df_fcg['_frags_u'].map(lambda x: frozenset(x))

    # extract edge informations for references
    df_fcg_ref["edges"] = df_fcg_ref["_fcg"].map(lambda x: x.edges(data=True))
    df_fcg_ref["edges"] = df_fcg_ref["edges"].map(lambda x: filter_edges_attributes(x, data))
    df_fcg_ref['_frags_u'] = df_fcg_ref['_frags_u'].map(lambda x: frozenset(x))

    # run annotation
    df_fcg['_pnp_ref'] = df_fcg.apply(lambda x: get_pnp_references(x['edges'], df_fcg_ref, x['_frags_u']), axis=1)
    df_fcg['pnp_fcg'] = df_fcg['_pnp_ref'].map(lambda x: True if len(x) == 0 else False)
    df_fcg['pnp_mol'] = df_fcg.groupby('idm')[['pnp_fcg']].transform(lambda x: True if any(x) else False)

    return df_fcg
