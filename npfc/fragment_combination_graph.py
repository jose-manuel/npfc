"""
Module fragment_combination_graph
=================================
This modules contains the functions for generating fragments combination graphs
from individual molecules:

It provides functions to annotate pseudo-Natural Products as well.
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
# graph
import networkx as nx
# docs
from typing import List
# dev
from npfc import fragment_combination_point
from npfc import draw
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# for computing fragment graphs
DF_FG_COLS = ['idm',
              'inchikey',
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
              '_comb',
              '_comb_u',
              'fcg_str',
              '_d_aidxs',
              '_colormap',
              '_fcg',
              'mol',
              '_d_mol_frags',
              '_d_fcp_labels',
              ]

# for annotating fragment graphs with PNP
DF_PNP_COLS = DF_FG_COLS + ['pnp_fm', 'pnp_mol']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _clear_ffs(df_fcc: DataFrame) -> DataFrame:
    """Clear ffs combinations by discarding any combination in which the smaller
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
        G = nx.from_pandas_edgelist(df_fcc_clean, "fid1", "fid2")  # we do not care for parallel edges here, just for the node ids in each subgraphs
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

    use_fcp_labels = True
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
            d_fcp_labels = {}
            for j in range(len(df_fcc_clean.index)):
                row = df_fcc_clean.iloc[j]
                # idf1
                if row["idf1"] not in d_aidxs.keys():
                    d_aidxs[row["idf1"]] = [row["_aidxf1"]]
                    d_frags[row["idf1"]] = row["mol_frag_1"]
                    try:
                        d_fcp_labels[row['idf1']] = row['_fcp_labels_1']
                    except KeyError:
                        use_fcp_labels = False
                elif row["_aidxf1"] not in d_aidxs[row["idf1"]]:
                    d_aidxs[row["idf1"]].append(row["_aidxf1"])
                    try:
                        d_fcp_labels[row['idf1']] = row['_fcp_labels_1']
                    except KeyError:
                        use_fcp_labels = False
                # idf2
                if row["idf2"] not in d_aidxs.keys():
                    d_aidxs[row["idf2"]] = [row["_aidxf2"]]
                    d_frags[row["idf2"]] = row["mol_frag_2"]
                    try:
                        d_fcp_labels[row['idf2']] = row['_fcp_labels_2']
                    except KeyError:
                        use_fcp_labels = False
                elif row["_aidxf2"] not in d_aidxs[row["idf2"]]:
                    d_aidxs[row["idf2"]].append(row["_aidxf2"])
                    try:
                        d_fcp_labels[row['idf2']] = row['_fcp_labels_2']
                    except KeyError:
                        use_fcp_labels = False

            # sort d_aidxs for reproducible colormaps  ### might be the cause of the wrong coloration in alternative fgraphs due to overlaps
            d_aidxs = OrderedDict(sorted(d_aidxs.items()))
            if use_fcp_labels:
                d_fcp_labels = OrderedDict(sorted(d_fcp_labels.items()))
            else:
                logging.warning('FCP labels are not available')

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
            # create two new columns: fcp_1 and fcp_2 for indicating fragment connection points of source and target
            df_fcc_clean['fcp_1'], df_fcc_clean['fcp_2'] = zip(*df_fcc_clean['fc'].map(lambda x: (x.split('[')[0].split('@')[1], x.split(']')[1].split('@')[1])))
            # compute the graph
            edge_attr = ['fcc', 'n_fcc', 'idm', 'fcp_1', 'fcp_2']
            edge_attr = [x for x in edge_attr if x in df_fcc_clean.columns]
            df_fcc_clean = df_fcc_clean.sort_values(['idf1', 'idf2'])
            df_fcc_edges = df_fcc_clean.copy()[['idf1', 'idf2'] + edge_attr]
            G = nx.from_pandas_edgelist(df_fcc_edges, source="idf1", target="idf2", edge_attr=edge_attr, create_using=nx.MultiGraph())  # multigraph instead of grpah to store parallel edges
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
            if not use_fcp_labels:
                d_fcp_labels = {}

            ds_fcg.append({'idm': gid, 'inchikey': g.iloc[0]['inchikey'], 'idfcg': str(i+1).zfill(3), 'nfrags': nfrags, 'nfrags_u': nfrags_u, 'ncomb': ncomb, 'ncomb_u': ncomb_u, 'hac_mol': hac_mol, 'hac_frags': hac_frags, 'perc_mol_cov_frags': perc_mol_cov_frags, '_frags': frags, '_frags_u': frags_u, '_comb': comb, '_comb_u': comb_u, 'fcg_str': fragment_combination_graph_str, '_d_aidxs': d_aidxs, '_colormap': colormap, '_fcg': G, 'mol': mol, '_d_mol_frags': d_frags, '_d_fcp_labels': d_fcp_labels})

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
    # hotfix: ####
    edges_ref = [tuple([set([er[0], er[1]]), er[2]]) for er in edges_ref]
    edges = [tuple([set([e[0], e[1]]), e[2]]) for e in edges]

    logging.debug('Reference edges: %s', edges_ref)
    for e in edges:
        logging.debug('current edge: %s', e)
        if e not in edges_ref:
            logging.debug('edge is not present within references!')
            return False
    return True


def filter_edges_attributes(edges: list, cols: list) -> list:
    """Networkx can either return one or all properties. Now that I am aware of this,
    I extract all attributes and then filter out the ones I do not want to use.

    :param edges: the edges as a list of tuple of syntax (u, v, d) with d being the dict with the attributes
    :param cols: the list of attributes to use for PNP labelling.
    """
    return [list([row[0], row[1], {k: v for k, v in row[2].items() if k in cols}]) for row in edges]


def get_pnp_references(edges: tuple, df_ref: DataFrame, target_node: frozenset = None) -> tuple:
    """Return a tuple of reference idms. A reference is recorded if all edges of the target molecule are included  at once
    within the reference row edges. In case no reference is found, an empty tuple is returned.

    :param edges: the edges as a list of tuple of syntax (u, v, d) with d being the dict with the attributes
    :param df_ref: the dataframe containing the edges to use for references.
    :param target_nodes: a frozenset of fragment ids found in the target edges to use for filering references to compare (optimization)
    """
    frags_u = frozenset([x[0] for x in edges] + [x[1] for x in edges])
    logging.debug('Input unique idfs to use: %s', frags_u)
    # avoid unnecessary graph comparisons: no way the syntehtic compound will be matching a NP fc graph if it does not even have all the nodes
    df_ref = df_ref[df_ref['_frags_u'].map(lambda x: frags_u.issubset(x))]
    logging.debug('Remaining reference entries after filtering by unique idfs: %s', len(df_ref))

    # get the subset of references that match with this particular molecule
    df_ref = df_ref[df_ref['edges'].map(lambda x: has_only_referenced_edges(edges, x))]

    # if no ref is found, just label as PNP
    if len(df_ref) < 1:
        return tuple()
    else:
        df_ref['idm_idfcg'] = df_ref['idm'].astype(str) + ':' + df_ref['idfcg'].astype(str)
        return tuple(df_ref['idm_idfcg'].values)


def annotate_pnp(df_fcg, df_fcg_ref, data=['fcc', 'fcp_1', 'fcp_2'], consider_symmetry=True) -> DataFrame:
    """Search and Identify for PNP molecules in the input DataFrame (df_fcg).
    PNP molecules are defined as molecules containing natural fragments combinations
    that are not found in a reference natural dataset.
    Fragment combinations are defined by extracting information from networkx graphs containing 3 informations:

        - source: id of fragment 1
        - target: id of fragment 2
        - attributes to consider: by default: fcc, fcp_1, fcp_2, etc.

    Three new columns are appended to the input DataFrame:

        - _pnp_ref: the list of a references found for the target fcg
        - pnp_fcg: True if the input fcg has no match with any reference fcg, False otherwise
        - pnp_mol: True if the input molecule has no matching fcg with any reference fcg, False otherwise

    When considering FCP, it is recommanded to consider symmetry. This results in ignoring the suffixes in FCPs,
    i.e. '1a' and '1b' become both '1'.

    :param df_fcg: the input DataFrame
    :param df_fcg_ref: the reference DataFrame
    :param data: the list of edge attributes to consider during fcg comparison
    :param consider_symmetry: consider fragment symmetry during annotating when using FCPs. To use FCPs data must include 'fcp_1' and/or 'fcp_2' (using or would not make sense here).
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

    # handle symmetry in case it is wished for
    if consider_symmetry and 'fcp_1' in data and 'fcp_2' in data:
        df_fcg['edges'] = df_fcg['edges'].map(fragment_combination_point.clear_fcp_suffixes_in_edges)
        df_fcg_ref['edges'] = df_fcg_ref['edges'].map(fragment_combination_point.clear_fcp_suffixes_in_edges)

    # run annotation
    df_fcg['_pnp_ref'] = df_fcg.apply(lambda x: get_pnp_references(x['edges'], df_fcg_ref, x['_frags_u']), axis=1)
    df_fcg['pnp_fcg'] = df_fcg['_pnp_ref'].map(lambda x: True if len(x) == 0 else False)
    df_fcg['pnp_mol'] = df_fcg.groupby('idm')[['pnp_fcg']].transform(lambda x: True if any(x) else False)

    # drop edges column since edges are included intact in G
    df_fcg = df_fcg.drop('edges', axis=1)

    return df_fcg


def regroup_edges_from_fcgs(df_fcg):
    """This function regroups all edges of a molecule. It can be applied to more than only one molecule at once.
    This will result in one single fcg per molecule, with all combinations and fragment occurrences counted only once,
    which is useful for reporting. (I should never have split fragment graphs).
    """
    idm_groups = df_fcg.groupby('idm')
    rows = []
    cols = ['idm', 'idfcg', 'idf1', 'fcp_1', 'fcc', 'idf2', 'fcp_2']
    for gid, g in idm_groups:
        combinations = set('-'.join(g['fcg_str']).split('-'))
        for c in combinations:
            f1 = c.split(':')[0]
            f2 = c.split(']')[1].split(':')[0]
            fcc = c.split('[')[1].split(']')[0]
            fcp1 = c.split('@')[1].split('[')[0]
            fcp2 = c.split('@')[-1]
            rows.append([gid, -1, f1, fcp1, fcc, f2, fcp2])
    df_edges = pd.DataFrame(rows, columns=cols)
    return df_edges


def get_ref_aidxs(df_fs):
    """Part of the hotfix for redundant FCGs.
    I did not record the occurrence id in the graphs, which was stupid.
    So now I need to use the df_fs to get the information instead.
    Needs to be used with fid col, which is defined in filter_out_fcgs_ffs_all.
    """
    return {k: v for k, v in zip(df_fs['fid'], df_fs['_aidxf'])}


def get_varying_d_aidxs(varying_fragments_occ, d):
    """Part of the hotfix for redundant FCGs.
    Attribute the atom indices to the fragmnet occurrences.
    """
    return {k: set(v) for k, v in d.items() if k in varying_fragments_occ}


def _keep_first_fcg(d1, d2):
    """Part of the hotfix for redundant FCGs.
    Determine if:

        0. both FCGs should be kept
        1. first FCG should be kept
        2. second FCG should be kept

    This is based on atom indices only.
    """
    c1 = 0
    c2 = 0
    # print(f"{d1=}")
    # print(f"{d2=}")
    for f1, a1 in d1.items():
        for f2, a2 in d2.items():
            # print(f"investigating: {f1} ({a1}) + {f2} ({a2})")
            if a1.issubset(a2):
                # print(f"{f1} included in {f2}")
                c1 += 1
            elif a2.issubset(a1):
                # print(f"{f2} included in {f1}")
                c2 += 1

    # print(f"{c1=}; {c2=}")

    if c1 == 0 and c2 > 0:
        return 1
    elif c2 == 0 and c1 > 0:
        return 2
    else:
        return 0


def filter_out_fcgs_ffs(df2, d):
    """Part of the hotfix for redundant FCGs.
    Apply on DF with FCGs of a same molecules.
    """
    df2 = df2.rename({'_frags': 'frags'}, axis=1)

    to_remove = []

    for i, row1 in enumerate(df2.itertuples()):
        frags1 = set(row1.frags)
        if row1.idfcg in to_remove:
            continue
        for row2 in df2.iloc[i+1:].itertuples():
            if row2.idfcg in to_remove:
                continue
            # print(f"comparing FCG {row1.idfcg} + FCG {row2.idfcg}")
            # print('\n' + "=" * 60 + '\n')
            # print(f"COMPARISON: FCG {row1.idfcg} + FCG {row2.idfcg}".center(60))
            # print('\n' + "=" * 60 + '\n')
            # define common / varying fragments
            frags2 = set(row2.frags)
            # print(f"{frags1=}")
            # print(f"{frags2=}")
            common_fragments_occ = frags1.intersection(frags2)
            varying_fragments_occ1 = frags1.difference(frags2)
            varying_fragments_occ2 = frags2.difference(frags1)

            # print("\ncommon_fragments_occ")
            # display(common_fragments_occ)
            # print("\nvarying_fragments_occ1")
            # display(varying_fragments_occ1)
            # print("\nvarying_fragments_occ2")
            # display(varying_fragments_occ2)

            # display(d)

            d1 = get_varying_d_aidxs(varying_fragments_occ1, d)
            d2 = get_varying_d_aidxs(varying_fragments_occ2, d)

            keep_first = _keep_first_fcg(d1, d2)
            if keep_first == 1:
                to_remove.append(row2.idfcg)
            elif keep_first == 2:
                to_remove.append(row1.idfcg)

            # print(f"keep_first: {keep_first}")
            # if keep_first == 1:
                # print(f"FCG {row2.idfcg} will be removed")
            # elif keep_first == 2:
                # print(f"FCG {row1.idfcg} will be removed")
            # else:
            #     pass
                # print("keep both FCG")

            # do not waste time on comparing other possibilities
            # when we already know that this FCG has to go
            if keep_first == 2:
                break

    # print('\n' + "=" * 60 + '\n')
    # print("RESULT".center(60))
    # print('\n' + "=" * 60 + '\n')

    to_remove = set(to_remove)
    # print("FCGs to remove:")
    # display(to_remove)

    return df2[~df2['idfcg'].isin(to_remove)]


def filter_out_fcgs_ffs_all(df_fcg, df_fs):
    """Part of the hotfix for redundant FCGs.
    """
    dfs = []
    df_fs['fid'] = df_fs['idf'].astype(str) + ':' + df_fs['idf_idx'].astype(str)
    for gid, g in df_fcg.groupby('idm'):
        d_ref = get_ref_aidxs(df_fs[df_fs['idm'] == gid])
        dfs.append(filter_out_fcgs_ffs(g, d_ref))

    return pd.concat(dfs)
