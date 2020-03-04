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
from copy import deepcopy
# data handling
import json
import base64
import numpy as np
from collections import OrderedDict
from itertools import chain
from collections import Counter
# chemoinformatics
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Mol
from rdkit.Chem import Atom
from rdkit.Chem import Bond
from rdkit.Chem import Draw
# 2D depiction of molecules
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdCoordGen
from pdbeccdutils.core.depictions import DepictionValidator
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
# dev
from npfc import utils


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
    In case the same color is found, a 15% darker shade is used after blending, so we can still distinguish common atoms and bonds.
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

        # in case of the same color being applied, use a 15% darker shade
        if new_color == old_color:
            new_color = tuple((x * 0.85 for x in color))
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


def clear_colors(mol: Mol, inplace: bool = True) -> Union[None, Mol]:
    """Clear any color highlight (atoms and bonds) on a molecule.
    This is equivalent to removing the _color and _num_colors atom and bond attributes.

    :param mol: the molecule to clear colors from
    :param inplace: if set to True, the molecule is directly modified. If set to False, a modified copy of the molecule is returned.
    :return: the modified molecule or nothing
    """
    props_to_clear = ['_color', '_num_colors']
    if not inplace:
        mol = deepcopy(mol)
    for p in props_to_clear:
        [a.ClearProp(p) for a in mol.GetAtoms()]
        [b.ClearProp(p) for b in mol.GetBonds()]
    if not inplace:
        return mol


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


def highlight_mol(mol: Mol,
                  colormap: 'ColorMap' = None,
                  output_file: str = None,
                  img_size: Tuple[int] = (300, 300),
                  debug: bool = False,
                  svg: bool = True,
                  legend: str = '') -> Image:
    """
    Draw an Image of a molecule with highlighted atoms and bonds according to a colormap.
    If no Colormap object is provided, no highlighting is done.

    .. image:: _images/draw_highlight.svg
        :align: center

    :param mol: the molecule to highlight
    :param colormap: the colormap to use for highlighting the molecule
    :param output_file: if speficied, the image is saved (format is deduced from extension)
    :param img_size: the size of the resulting Image
    :param debug: display atom indices on the structure
    :param svg: use SVG format instead of PNG
    :return: an Image of the highlighted molecule
    """
    if debug:
        mol = Mol(mol)
        [mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx())) for idx in range(mol.GetNumAtoms())]
    if colormap is None:
        highlightAtomList = []
        highlightAtomColors = {}
        highlightBondColors = {}
    else:
        highlightAtomList = [int(x) for x in list(colormap.atoms.keys())]
        highlightAtomColors = colormap.atoms
        highlightBondColors = colormap.bonds
    img = Draw.MolsToGridImage([mol],
                               molsPerRow=1,
                               subImgSize=img_size,
                               highlightAtomLists=[highlightAtomList],
                               highlightAtomColors=[highlightAtomColors],
                               highlightBondColors=[highlightBondColors],
                               useSVG=svg,
                               legends=[legend],
                               )
    # export img
    if output_file is not None:
        output_ext = output_file.split('.')[-1].upper()
        if output_ext == 'SVG' and not svg:
            raise ValueError(f"Error! output file extension is SVG but image format is PNG!")
        if output_ext == 'SVG':
            with open(output_file, 'w') as SVG:
                SVG.write(img.data)
        else:
            raise ValueError(f"Error! Unsupported extension '{output_ext}'!")
    return img


def highlight_mols(mols: List[Mol],
                   colormaps: List['ColorMap'] = [],
                   output_file: str = None,
                   sub_img_size: Tuple[int] = (300, 300),
                   max_mols_per_row: int = 5,
                   debug: bool = False,
                   svg: bool = True,
                   legends: List[str] = None,
                   ):
    """
    Draw an Image of a list of molecules with highlighted atoms and bonds
    according to a list of colormaps.

    :param mols: the molecules to highlight
    :param colormaps: the colormaps to use for highlighting the molecules
    :param output_file: if speficied, the image is saved (format is deduced from extension)
    :param sub_img_size: the size of the image of every molecule composing the grid
    :param max_mols_per_row: the maximum number of molecules displayed per row
    :param debug: display atom indices on the structure
    :param svg: use SVG format instead of PNG
    :return: an Image of the highlighted molecules
    """
    atom_lists = []
    colormaps_a = []
    colormaps_b = []

    if colormaps is None:
        for colormap in colormaps:
            atom_lists.append([])
            colormaps_a.append({})
            colormaps_b.append({})
    else:
        for colormap in colormaps:
            atom_lists.append([int(x) for x in list(colormap.atoms.keys())])
            colormaps_a.append(colormap.atoms)
            colormaps_b.append(colormap.bonds)

    if debug:
        mols = [Mol(mol) for mol in mols]
        for mol in mols:
            [mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx())) for idx in range(mol.GetNumAtoms())]

    img = Draw.MolsToGridImage(mols,
                               molsPerRow=max_mols_per_row,
                               subImgSize=sub_img_size,
                               highlightAtomLists=atom_lists,
                               highlightAtomColors=colormaps_a,
                               highlightBondColors=colormaps_b,
                               useSVG=svg,
                               legends=legends,
                               )
    # export img
    if output_file is not None:
        output_ext = output_file.split('.')[-1].upper()
        if output_ext == 'SVG' and not svg:
            raise ValueError(f"Error! output file extension is SVG but image format is PNG!")
        if output_ext == 'SVG':
            with open(output_file, 'w') as SVG:
                SVG.write(img.data)
        else:
            raise ValueError(f"Error! Unsupported extension '{output_ext}'!")
    return img


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


def fc_graph(fc_graph: Graph, colormap_nodes: List[Tuple[float]] = None, output_file: str = None) -> Figure:
    """
    Return a matplotlib Figure of a networkx graph.

    :param fc_graph: a networkx Graph object of the fragment combinations
    :param colormap_nodes: a colormap of RGB values for the nodes (i.e. [(0, 0, 1), (0, 1, 0)])
    :return: a matplotlib Figure object
    """
    if isinstance(fc_graph, base64.bytes_types):
        fc_graph = utils.decode_object(fc_graph)

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
            alpha=1,
            font_size=16,
            with_labels=True,
            connectionstyle='arc3,rad=0.9',
            )
    nx.draw_networkx_edge_labels(fc_graph,
                                 pos,
                                 edge_labels=edges_info,
                                 font_color='red',
                                 font_size=14,
                                 )
    if output_file is not None:
        output_file_format = output_file.split('.')[-1].upper()
        plt.savefig(output_file, format=output_file_format)
    plt.close()
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
        colormap_nodes = list(colormap.fragments.values())

    return fc_graph(row["fc_graph"], colormap_nodes=colormap_nodes)


def rescale(mol: Mol, f: float = 1.4):
    """Rescale the coordinates of a Mol with a factor f.

    :param mol: a Mol which is modified in place.
    :param f: the factor for rescaling coordinates
    """
    tm = np.zeros((4, 4), np.double)
    for i in range(3):
        tm[i, i] = f
    tm[3, 3] = 1.0
    Chem.TransformMol(mol, tm)


def depict_mol(mol: Mol, methods: List[str] = ["CoordGen", "rdDepictor"], consider_input: bool = True) -> Mol:
    """
    Returns the "best" 2D depiction of a molecule according the methods in METHODS_2D.
    Currently two methods are available:

        - CoordGen
        - rdDepictor

    A perfect score of 0 means the depiction is good enough (no overalapping atom/bonds)
    and it is not worth computing other depictions. When no perfect score is reached,
    the depiction with lowest score is retrieved. In case of tie, the first method applied
    is preferred. In case the input molecule contains input coordinates, they can be compared
    to the methods as 'Input' (lowest priority).

    The method used for depicting the molecule is stored as molecule property: "_2D".

    For CoordGen, 2D representations are automatically rescaled with a factor of 1.4.

    Examples of 2D coordinates computed with listed methods (source: SDF from ChEMBL)

    - Simple case

    .. image:: _images/draw_simple.png
        :align: center

    Most of molecules including macrocycles are usually better rendered with CoordGen.

    - Medium case

    .. image:: _images/draw_medium.png
        :align: center

    In same cases, Avalon performs better than CoordGen.

    - Complex Case

    .. image:: _images/draw_hard.png
        :align: center

    .. note:: Avalon was however removed from the available methods as it produces mostly errors since latest RDKit update...

    For some molecules, none of the methods yield a "perfect score". The depiction with the lowest score is thus selected.

    :param mol: the input molecule
    :param methods: a list of methods to apply. Currently supported: CoordGen, rdDepictor.
    :param consider_input: consider the input coordinates (if any), for determining the best 2D representation
    :return: the molecule with 2D coordinates and a new "_2D" property with the information of which depictor was selected.
    """

    # methods
    METHODS = {'CoordGen': lambda x: rdCoordGen.AddCoords(x),
               'rdDepictor': lambda x: rdDepictor.Compute2DCoords(x),
               'Input': lambda x: x,
               }

    depictions = OrderedDict()

    for method in methods:
        # copy the input mol so input coordinates are not modified
        depiction_mol = Mol(mol)
        # compute the depiction (in place)
        METHODS[method](depiction_mol)
        # coordgen creates very small 2D representations, so let's rescale them
        if method == "CoordGen":
            rescale(depiction_mol)
        # score the depiction
        dv = DepictionValidator(depiction_mol)
        [a.SetProp('name', str(a.GetIdx())) for a in dv.mol.GetAtoms()]  # bug fix for version 5+,add name property to atoms
        depiction_score = dv.depiction_score()

        # exit if perfect score, record depiction for selection otherwise
        if depiction_score == 0:
            depiction_mol.SetProp("_2D", method)
            return depiction_mol
        depictions[method] = (depiction_score, depiction_mol)

    # no perfect score was reached until now, so test input coordinates if any
    if consider_input and mol.GetNumConformers() > 0:
        method = "Input"
        dv = DepictionValidator(mol)
        [a.SetProp('name', str(a.GetIdx())) for a in dv.mol.GetAtoms()]  # bug fix for version 5+,add name property to atoms
        depiction_score = dv.depiction_score()
        if depiction_score == 0:
            mol.SetProp("_2D", method)
            return mol
        depictions[method] = (depiction_score, mol)
    elif consider_input and mol.GetNumConformers() == 0:
        logging.debug("No input coordinates to use for Input method, so skipping it!")

    # retrieve best depiction possible
    best_method = min(depictions, key=lambda k: depictions[k][0])
    best_depiction_mol = depictions[best_method][1]
    best_depiction_mol.SetProp("_2D", best_method)

    return best_depiction_mol


def draw_reaction(mol1: Mol, mol2: Mol, sub_img_size=(200, 200), svg=True):
    """Wrapper function around RDKit ReactionToImage function.

    :param mol1: the molecule to display left
    :param mol2: the molecule to display right
    :param sub_img_size: the size of the molecules
    :param svg: return the image in SVG text, return a PIL Image otherwise.
    :return: the reaction as an image.
    """
    if isinstance(mol1, Mol):
        mol1 = Chem.MolToSmiles(mol1)
    if isinstance(mol2, Mol):
        mol2 = Chem.MolToSmiles(mol2)
    rxn_str = f"{mol1}>>{mol2}"
    logging.debug(f"rxn_str={rxn_str}")
    rxn = rdChemReactions.ReactionFromSmarts(rxn_str, useSmiles=True)
    img = Draw.ReactionToImage(rxn, subImgSize=sub_img_size, useSVG=svg)

    return img


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ColorMap:
    """A class containing all the required information for highlighting a Molecule.
    It is represented by the count of fragments, atom- and bond colors.
    """

    def __init__(self, mol: Mol, d_aidxs: OrderedDict, colors: OrderedDict = colors, reset_colors: bool = False):
        """
        :param mol: the molecule to highlight. Atom/Bond properties '_color' and 'num_colors' are modified in place.
        :param d_aidxs: a dictionary containing fragment ids as keys and molecule atom indices as values ({fid: [aidxs]})
        :param colors: the color palette to use, by default the draw.colors palette is used (red, green, blue, orange and purple).
        :param reset_colors: clear any pre-existing colors in atoms and bonds. This ensures independent colorations using a same molecule.
        """
        self.fragments, self.atoms, self.bonds = self._compute_colormap(mol, d_aidxs, colors, reset_colors)

    def __repr__(self):
        """Return a string representation of the ColorMap object
        """
        # listing a huge dictionary of tuples is not great for representing data,
        # so I just list the number of different colors found
        num_frags = len(list(self.fragments.keys()))
        num_atom_colors = len(set(self.atoms.values()))
        # bonds require special handling because of a hack
        bond_colors = set(self.bonds.values())
        try:
            bond_colors.remove((1, 1, 1))  # do not count hard-coded white bonds, also white color should never happen during blending anyway
        except KeyError:
            pass
        num_bond_colors = len(bond_colors)

        return str(f"ColorMap(nf={num_frags}, nac={num_atom_colors}, nbc={num_bond_colors})")

    def _compute_colormap(self, mol: Mol, d_aidxs: Dict, colors, reset_colors):
        """
        Compute a colormap for highlighting a molecule given a dictionary of fragments {fid: [aidxs]} and specified RGB colors.

        The following rules are applied for consistant highilighting:

                 - only one color is attributed per fragment type (defined by fragment id, i.e. f1 in f1:0)
                 - in case there are more fragments than colors, a same color can be used for several fragments
                 - when 2 fragments of different colors overlap, their colors are blended on overlapping atoms/bonds
                 - when 2 fragments of same colors overlap, a 15% darker shade is used on overlapping atoms/bonds, so these can be distinguished

        :param mol: the molecule to highlight
        :param d_aidxs: a dictionary of fragments atom indices with fragment ids as keys and list of iterable as values
        :param colors: a color palette
        :param reset_colors: reset colors found on a molecule (inplace only)
        """
        if reset_colors:
            clear_colors(mol, inplace=True)
        colors_k = list(colors.keys())  # colors is a OrderedDict, so no need for resorting colors
        aidxs_colored = set()  # atoms
        bidxs_colored = set()  # bonds
        colormap_f = OrderedDict()  # fragments are colored and stored in order
        # for i, (fid, aidxs_l) in enumerate(zip(fids, l_aidxs)):
        for i, (fid, aidxs_l) in enumerate(d_aidxs.items()):
            # pick a color
            color = colors[colors_k[i % len(colors_k)]]
            colormap_f[fid] = color
            # print(f"fid: {fid}, color: {color}, aidxs_l: {aidxs_l}")
            # highlight the corresponding fragment
            for aidxs in aidxs_l:
                set_atoms_color(mol, aidxs, color)
                bidxs = get_bidxs(mol, aidxs)
                set_bonds_color(mol, bidxs, color)
                aidxs_colored.update(set(aidxs))
                bidxs_colored.update(set(bidxs))

        colormap_a = get_atoms_color(mol, aidxs_colored)  # atoms
        colormap_b = get_bonds_color(mol, bidxs_colored)  # bonds

        # highlight bonds not colored as white, since we get sometimes the RDKit default highlight otherwise
        # this might not be the best idea ever, but I could not find any hidden property on the mol, atoms or bonds
        # responsible for this default for behavior. To compensate, I tried to make it as fast as possible,
        # hence the cryptic writing below.
        # Also I do not update the _color and _num_colors properties, so further color blending should not be an issue
        bidxs_white = {bidx: (1, 1, 1) for bidx in set([b.GetIdx() for b in mol.GetBonds()]) - set(list(colormap_b.keys()))}
        #                      white                            all bidx of the mol          -        all colored bidx
        # update colormap with white bond indices
        logging.debug('Highlighting out-of-fragment-bonds in white')
        colormap_b.update(bidxs_white)

        # return results as tuple
        return (colormap_f, colormap_a, colormap_b)
