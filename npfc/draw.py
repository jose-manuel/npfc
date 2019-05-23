"""
Module draw
==============
This modules contains methods for drawing images of molecules with
highlighted fragments.
A special care was given to blending colors for overlapping fragments.
"""

# standard
import logging
from math import sqrt
# data handling
import json
import base64
import pickle
from collections import OrderedDict
from itertools import chain
from collections import Counter
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import Atom
from rdkit.Chem import Bond
from rdkit.Chem import Draw
# graph
import networkx as nx
# docs
import matplotlib.pyplot as plt  # required for creating a canvas for displaying graphs
from matplotlib.figure import Figure
from networkx.classes.graph import Graph
from PIL.Image import Image
from typing import Union
from typing import Set
from typing import List
from typing import Tuple
from typing import Dict
from pandas import Series


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# OrderedDict so fragment 1 is always colored in red, fragment 2 in green, etc.
colors = OrderedDict()
colors["red"] = (1.0, 0.6, 0.6)
colors["green"] = (0.2, 1.0, 0.2)
colors["blue"] = (0.4, 0.6, 1.0)
colors["orange"] = (0.9569, 0.6667, 0.2588)
colors["purple"] = (0.8392, 0.6275, 0.7686)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _blend_color_value(val1: int, val2: int, alpha: float) -> Tuple[int]:
    return sqrt((1-alpha) * (val1 ** 2) + alpha * (val2 ** 2))


def blend_color(color1: Tuple[int], color2: Tuple[int], alpha: float = 0.5) -> Tuple[int]:
    """Blend two colors after Downgoat'reply at:

    https://stackoverflow.com/questions/726549/algorithm-for-additive-color-mixing-for-rgb-values

    The formala below is applied to each color channel (R, G, and B):

    .. math::
        blended = \\sqrt{(1 - alpha) * val1^2 + alpha * val2^2}

    with val1: color channel in color1; val2: corresponding color channel in color2; alpha: transparency factor

    :param color1: the first color to blend
    :param color2: the second color to blend
    :param alpha: the blending factor, 0.5 means as much as color1 than color2
    """
    return tuple((_blend_color_value(x1, x2, alpha) for x1, x2 in zip(color1, color2)))


def get_bidxs(mol: Mol, aidxs: Set[int]) -> Set[int]:
    """From an iterable of atom indices (aidxs), return a set of bond indices (bidxs)
    found between atoms in aidxs.

    :param mol: the molecule to highlight
    :param aidxs: a set with the atom indices of the fragment to highlight
    :return: a set with the corresponding bond indices
    """
    bidxs = []
    for aidx in aidxs:
        bonds = [b for b in mol.GetAtomWithIdx(aidx).GetBonds()]
        bidx = [b.GetIdx() for b in bonds]
        bidxs.append(bidx)
    # flatten the list of bond indices (bidx)
    bidxs_merged = list(chain.from_iterable(bidxs))  # all bonds
    # count number of times each bidx is found
    counter = Counter(bidxs_merged)
    # keep only bidx that appear > 1
    return set({x: counter[x] for x in counter if counter[x] > 1}.keys())


def get_bidx_all(mol: Mol, aidxs: Tuple[int]) -> Set[int]:
    """From an iterable of atom indices (aidxs), return a set of all bond indices (bidxs)
    connected to a least one specified atom.

    :param mol: the molecule to highlight
    :param aidxs: a set with the atom indices of the fragment to highlight
    :return: a set with the corresponding bond indices
    """
    bidx = []
    for aidx in aidxs:
        [bidx.append(b.GetIdx()) for b in mol.GetAtomWithIdx(aidx).GetBonds()]
    return set(bidx)


def set_atom_or_bond_color(atom_or_bond: Union[Atom, Bond], color: Tuple[float]):
    """Set a color (RGB , i.e. (0, 0, 1) as str property using JSON to an atom or a bond of a RDKit Molecule.
    This is possible because they share the same API for doing executing this task.
    If a color is already defined, then an attempt is made to blend it with
    the new color. The number of colors that were already mixed is used for weighting the blending.
    In case the same color is found, a 10% darker shade is used after blending, so we can still distinguish common atoms and bonds.
    The atom or bond property is modified in place (_num_colors, _color).

    .. note:: Due to the color palette I am using, I could only get a "olive green" when blending red and green. To get a proper golden yellow, I hard-coded the color for this particular blending.

    :param atom_or_bond: either an atom or a bond
    :param color: a tuple of RGB values
    """
    if atom_or_bond.HasProp("_num_colors"):
        num_colors = int(atom_or_bond.GetProp("_num_colors"))
    else:
        num_colors = 0

    # if already a color on the atom, blend it with the new one
    if num_colors > 0:
        old_color = get_atom_or_bond_color(atom_or_bond)
        # blend with the old color only if it is not the same color
        old_alpha = 1 - (1/(num_colors + 1))
        alpha = 1 - old_alpha
        new_color = blend_color(old_color, color, alpha)
        # when I mix red and green, I want to get some yellow, not some olive color,
        # so here I hardcode the yellow I long for (golden yellow)
        if tuple(round(x, 4) for x in new_color) == (0.7211, 0.8246, 0.4472):
            new_color = (1.0, 0.9294, 0.0)

        if new_color == old_color:
            new_color = tuple((x * 0.9 for x in color))
    else:
        new_color = color

    atom_or_bond.SetProp("_color", json.dumps(list(new_color)))
    atom_or_bond.SetProp("_num_colors", str(num_colors+1))


def get_atom_or_bond_color(atom_or_bond: Union[Atom, Bond]) -> Tuple[float]:
    if atom_or_bond.HasProp("_color"):
        return tuple(json.loads(atom_or_bond.GetProp("_color")))


def set_atoms_color(mol: Mol, aidxs: Set[int], color: Tuple[float]):
    """Set all specified atoms in a given color. If some atoms are already set to
    another color, then an attempt is made to blend colors.

    :param mol: the molecule to highlight
    :param aidxs: the indices of the atoms to highlight
    :param color: a color as a tuple of R, G, B values within a range of [0;1]
    """
    atoms = [mol.GetAtomWithIdx(x) for x in aidxs]
    [set_atom_or_bond_color(a, color) for a in atoms]


def get_atoms_color(mol: Mol, aidxs: Set[int]) -> Dict:
    atoms = [mol.GetAtomWithIdx(x) for x in aidxs]
    return {a.GetIdx(): get_atom_or_bond_color(a) for a in atoms}


def set_bonds_color(mol: Mol, bidxs: Set[int], color: Tuple[float]):
    """Set all specified bonds in a given color. If some bonds are already set to
    another color, then an attempt is made to blend colors.

    :param mol: the molecule to highlight
    :param bidxs: the indices of the bonds to highlight
    :param color: a color as a tuple of R, G, B values within a range of [0;1]
    """
    bonds = [mol.GetBondWithIdx(x) for x in bidxs]
    [set_atom_or_bond_color(a, color) for a in bonds]


def get_bonds_color(mol: Mol, bidxs: Set[int]) -> Dict:
    bonds = [mol.GetBondWithIdx(x) for x in bidxs]
    return {b.GetIdx(): get_atom_or_bond_color(b) for b in bonds}


def scale_rgb(rgb_color: Tuple[int]) -> Tuple[float]:
    """Scale down a RGB value (0, 0, 255) to (0.0, 0.0, 1.0).
    Useful for trying new colors from the web as RDKit only accepts [0;1] range.

    :param rgb_color: a color with RGB values in range [0;250]
    :return: a color with RGB values in range [0;1]
    """
    if any((isinstance(v, int) or v > 1.0) for v in rgb_color):
        return tuple([round(x/255, 4) for x in rgb_color])


def scale_rgb_colormap(colormap: Dict) -> Dict:
    """Scale down the RGB values in a colormap so that they are in acceptable range
    for RDKit (i.e. [0;1]).

    :param colormap: a dictionary of syntax: atom_or_bond_idx: (R, G, B)
    :return: the same colormap but scaled down.
    """
    return {k: scale_rgb(colormap[k]) for k in colormap.keys()}


def _compute_colormap(mol: Mol, d_aidxs: Dict, colors):
    """
    Compute a colormap for highlighting a molecule given a dictionary of fragments {fid: [aidxs]} and specified RGB colors.

    The following rules are applied for consistant highilighting:

             - only one color is attributed per fragment
             - in case there are more fragments than colors, a same color can be used for several fragments
             - when 2 fragments of different colors overlap, their colors are blended on overlapping atoms/bonds
             - when 2 fragments of same colors overlap, a 10% darker shade is used on overlapping atoms/bonds, so these can be distinguished

    :param mol: the molecule to highlight
    :param d_aidxs: a dictionary of fragments atom indices
    :param colors: a color palette
    """
    # extract all fids and sort them by alphabetical order for reproducible coloring
    fids = list(d_aidxs.keys())
    fids.sort()
    # get corresponding values
    l_aidxs = [d_aidxs[fid] for fid in fids]
    colors_k = list(colors.keys())  # colors is a OrderedDict, so no need for resorting colors
    colormap_a = {}  # atoms
    colormap_b = {}  # bonds
    colormap_f = OrderedDict()  # fragments are colored and stored in order
    for i, (fid, aidxs) in enumerate(zip(fids, l_aidxs)):
        # pick a color
        color = colors[colors_k[i % len(colors_k)]]
        colormap_f[fid] = color
        # highlight the corresponding fragment
        set_atoms_color(mol, aidxs, color)
        bidxs = get_bidxs(mol, aidxs)
        set_bonds_color(mol, bidxs, color)
        # colormaps for the current fragment, might not be the same as input if atom/bond was already colored
        colormap_a_curr = get_atoms_color(mol, aidxs)
        colormap_b_curr = get_bonds_color(mol, bidxs)
        # colormaps for the whole molecule
        colormap_a.update(colormap_a_curr)
        colormap_b.update(colormap_b_curr)
        # display detailed results in case of debugging
        if logging.getLogger().level == logging.DEBUG:
            padding = 30
            logging.debug(f"Highlighting the molecule #{i}:")
            print("=" * padding + "\n" + "highlightAtomLists".center(padding) + "\n" + "=" * padding)
            print([a for a in list(colormap_a.keys())])
            print("=" * padding + "\n" + "highlightAtomColors".center(padding) + "\n" + "=" * padding)
            [print(f"{k}: {v}") for k, v in colormap_a.items()]
            print("=" * padding + "\n" + "highlightBondColors".center(padding) + "\n" + "=" * padding)
            [print(f"{k}: {v}") for k, v in colormap_b.items()]

    # highlight bonds not colored as white, since we get sometimes the RDKit default highlight otherwise
    # this might not be the best idea ever, but I could not find any hidden property on the mol, atoms or bonds
    # responsible for this default for behavior. To compensate, I tried to make it as fast as possible,
    # hence the cryptic writing below.
    # Also I do not update the _color and _num_colors properties, so further color blending should not be an issue
    bidxs_white = {bidx: (1, 1, 1) for bidx in set([b.GetIdx() for b in mol.GetBonds()]) - set(list(colormap_b.keys()))}
    #                      white                            all bidx of the mol          -        all colored bidx
    # update colormap with white bond indices
    logging.debug(f"")
    colormap_b.update(bidxs_white)

    # return results as dict
    return {'atoms': colormap_a, 'bonds': colormap_b, 'fragments': colormap_f}


# def draw_mols_frags(mols: List[Mol],
#                     ll_aidxs: List[List[Set[int]]],
#                     colors: List[Tuple[float]] = colors,
#                     debug: bool = False,
#                     size: Tuple[int] = (300, 300),
#                     max_mols_per_row: int = 5,
#                     ) -> Image:
#     """Draw a PNG grid image of molecules with highlighted fragments using the RDKit Drawing functionalities.
#     For each molecule, several fragments can be highlighted. In case a same atom or bond is highlighted several times,
#     an attempt is made to blend colors.
#     A debug mode can be activated for displaying atom indices on the image, in addition to the highlighting.
#     If no fragment indices are specified, no highlighting is done.
#
#     The following rules are applied for consistant highilighting:
#
#         - only one color is attributed per fragment
#         - in case there are more fragments than colors, a same color can be used for several fragments
#         - when 2 fragments of different colors overlap, their colors are blended on overlapping atoms/bonds
#         - when 2 fragments of same colors overlap, a 10% darker shade is used on overlapping atoms/bonds, so these can be distinguished
#
#     .. note:: The parameter name ll_aidxs might seem confusing. Let's say we start with a fragment of 3 atoms, one atom index (aidx) might be for instance 0, 1 or 2. The whole fragment is thus represented with aidxs={0, 1, 2}. The list of all fragments for a molecule is then l_aidxs=[{0, 1, 2}]. And if there are several molecules, then we need several lists of l_aidxs as well, hence ll_aidxs.
#
#     :param mols: the molecules to draw
#     :param ll_aidxs: a list of list of sets of atom indices. Each set corresponds to a fragment to highlight.
#     :param debug: display atom indices on the image too
#     :param size: the image size
#     :param max_mols_per_row:
#     :return: a PNG image of the grid of molecules with highlighted fragments
#     """
#     if len(mols) != len(ll_aidxs):
#         raise ValueError(f"Error! The number of mols is different than the number of fragment lists ({len(mols)} != {len(ll_aidxs)})")
#
#     atom_lists = []
#     colormaps_a = []
#     colormaps_b = []
#
#     for i, mol in enumerate(mols):
#         if debug:
#             mol = Mol(mol)
#             [mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx())) for idx in range(mol.GetNumAtoms())]
#
#         l_aidxs = ll_aidxs[i]
#         colormap = _compute_colormap(mol, l_aidxs, colors)
#         atom_lists.append(list(colormap['atoms'].keys()))
#         colormaps_a.append(colormap['atoms'])
#         colormaps_b.append(colormap['bonds'])
#
#     return Draw.MolsToGridImage(mols,
#                                 molsPerRow=max_mols_per_row,
#                                 subImgSize=size,
#                                 highlightAtomLists=atom_lists,
#                                 highlightAtomColors=colormaps_a,
#                                 highlightBondColors=colormaps_b,
#                                 )


def highlight_mol(mol: Mol, colormap: Dict, img_size: Tuple[int] = (300, 300), debug: bool = False) -> Image:
    """
    Draw an Image of a molecule with highlighted atoms and bonds according to a colormap.

    :param mol: the molecule to draw
    :param colormap: a dictionary containing 3 keys: 'fragments', 'atoms' and 'bonds' and associated colors
    :param img_size: the size of the resulting Image
    :param debug: display atom indices on the structure
    :return: a PNG Image of the highlighted molecule
    """
    if debug:
        mol = Mol(mol)
        [mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx())) for idx in range(mol.GetNumAtoms())]

    return Draw.MolsToGridImage([mol],
                                molsPerRow=1,
                                subImgSize=img_size,
                                highlightAtomLists=[[int(x) for x in list(colormap['atoms'].keys())]],
                                highlightAtomColors=[colormap['atoms']],
                                highlightBondColors=[colormap['bonds']],
                                )


def highlight_mols(mols, colormaps, sub_img_size=(300, 300), max_mols_per_row=5, debug=False):
    atom_lists = []
    colormaps_a = []
    colormaps_b = []

    for colormap in colormaps:
        atom_lists.append([int(x) for x in list(colormap['atoms'].keys())])
        colormaps_a.append(colormap['atoms'])
        colormaps_b.append(colormap['bonds'])

    if debug:
        mols = [Mol(mol) for mol in mols]
        for mol in mols:
            [mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx())) for idx in range(mol.GetNumAtoms())]

    return Draw.MolsToGridImage(mols,
                                molsPerRow=max_mols_per_row,
                                subImgSize=sub_img_size,
                                highlightAtomLists=atom_lists,
                                highlightAtomColors=colormaps_a,
                                highlightBondColors=colormaps_b,
                                )


def _get_edge_info(fc_graph: Graph) -> Dict:
    """
    Use the first associated data of edges for edge labelling of a networkx graph.

    :param fc_graph: a Fragment Combination graph
    :return: a Dict of syntax {(node1, node2): data}
    """
    d = {}
    for edge in list(fc_graph.edges(data=True)):
        d[(edge[0], edge[1])] = list(edge[2].values())[0]
    return d


def fc_graph(fc_graph: Graph, colormap_nodes: List[Tuple[float]] = None) -> Figure:
    """
    Return a matplotlib Figure of a networkx graph.

    :param fc_graph: a networkx Graph object of the fragment combinations
    :param colormap_nodes: a colormap of RGB values for the nodes (i.e. [(0, 0, 1), (0, 1, 0)])
    :return: a matplotlib Figure object
    """
    if isinstance(fc_graph, base64.bytes_types):
        fc_graph = pickle.loads(base64.b64decode(fc_graph))

    if colormap_nodes is None:
        # define a 2D list instead of a single tuple to avoid matplotlib warning
        colormap_nodes = [(0.7, 0.7, 0.7)] * len(list(fc_graph.nodes()))

    pos = nx.spring_layout(fc_graph)
    edges_info = _get_edge_info(fc_graph)
    figure = plt.figure()
    nx.draw(fc_graph,
            pos,
            edge_color='black',
            width=1,
            linewidths=1,
            node_size=2000,
            node_color=colormap_nodes,
            alpha=0.95,
            with_labels=True,
            )
    nx.draw_networkx_edge_labels(fc_graph,
                                 pos,
                                 edge_labels=edges_info,
                                 font_color='red',
                                 )
    return figure


def fc_graph_from_series(row: Series, colormap_nodes_name: str = None) -> Figure:
    """
    Return a matplotlib Figure of a networkx graph.

    :param row: a row from a Fragment Map DataFrame (df_map)
    :param colormap_nodes_name: the name of the Series axe (column) from which to retrieve the nodes colormap
    :return: a matplotlib Figure object
    """
    if colormap_nodes_name is None:
        colormap_nodes = None
    elif colormap_nodes_name == "fid":
        colormap = row["colormap"]
        colormap_nodes = list(colormap["fragments"].values())

    return fc_graph(row["fc_graph"], colormap_nodes=colormap_nodes)
