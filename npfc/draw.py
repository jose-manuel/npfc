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
import math
# data handling
import json
import base64
import numpy as np
from collections import OrderedDict
from itertools import chain
from collections import Counter
# chemoinformatics
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Mol
from rdkit.Chem import Atom
from rdkit.Chem import Bond
from rdkit.Chem import Draw
# 2D depiction of molecules
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdCoordGen
# graph
import networkx as nx
# docs
import matplotlib
import matplotlib.pyplot as plt  # required for creating a canvas for displaying graphs
from matplotlib.figure import Figure
from networkx.classes.graph import Graph
import seaborn as sns
from PIL.Image import Image
from typing import Union
from typing import Set
from typing import List
from typing import Tuple
from typing import Dict
# dev
from npfc import utils
from npfc import fragment_combination_graph
# tmp
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdCoordGen
from scipy.spatial import KDTree


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# in python 3.7+ dict are odered so red will always #1, green #2, etc.
# but I got a bug so now I keep using a nice defined OrderedDict
DEFAULT_PALETTE = OrderedDict()
DEFAULT_PALETTE['red'] = (1.0, 0.6, 0.6)
DEFAULT_PALETTE['green'] = (0.2, 1.0, 0.2)
DEFAULT_PALETTE['blue'] = (0.4, 0.6, 1.0)
DEFAULT_PALETTE['orange'] = (0.9569, 0.6667, 0.2588)
DEFAULT_PALETTE['purple'] = (0.8392, 0.6275, 0.7686)
DEFAULT_PALETTE['yellow'] = (1.0, 0.9294, 0.0)
DEFAULT_PALETTE['teal'] = (0.5725, 0.9608, 0.9882)
DEFAULT_PALETTE['gray'] = (0.7294, 0.7294, 0.7294)
#
#
#
# DEFAULT_PALETTE = {'red': (1.0, 0.6, 0.6),
#                    'green': (0.2, 1.0, 0.2),
#                    'blue': (0.4, 0.6, 1.0),
#                    'orange': (0.9569, 0.6667, 0.2588),
#                    'purple': (0.8392, 0.6275, 0.7686),
#                    'yellow': (1.0, 0.9294, 0.0),
#                    'teal': (0.5725, 0.9608, 0.9882),
#                    'gray': (0.7294, 0.7294, 0.7294),
#                    }

# matplotlib.colors.to_rgb('#FF0000')
# matplotlib.colors.to_hex(((1.0, 0.0, 0.0)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def get_d_aidxs_for_rings(mol: Mol, fuse_rings: bool = False) -> dict:
    """Return a dictionary of atom indices (d_aidxs) for highlighting rings in molecules.

    This is purely a helper function to instanciate a ColorMap object, which can highlight SSSR or ring systems instead of fragments.

    :param mol: the input molecule
    :param fuse_rings: if False: highlight SSSR, if True: highlight fused ring systems
    :return: a dictionary of syntax: {"R#0": [(0, 1, 2, 3, 4)], "R#1": [(5, 6, 7, 8, 9)]}. It is a list of tuple for compatibility reasons with fragment highlighting. It could also be edited to highlight rings of a certain size, etc.
    """
    ring_atoms_l = mol.GetRingInfo().AtomRings()
    if fuse_rings:
        ring_atoms_l = utils.fuse_rings(ring_atoms_l)
    return {f"R#{i}": [ring_atoms] for i, ring_atoms in enumerate(ring_atoms_l)}


def mol(mol: Mol,
        colormap: 'ColorMap' = None,
        output_file: str = None,
        img_size: Tuple[int] = (400, 400),
        debug: bool = False,
        svg: bool = True,
        legend: str = '') -> Image:
    """
    Draw an Image of a molecule with highlighted atoms and bonds according to a colormap.
    If no Colormap object is provided, no highlighting is done.

    This code is based on the 2020.03 RDKit release and the picture below is not yet updated.

    .. image:: _images/draw_highlight.svg
        :align: center

    :param mol: the molecule to highlight
    :param colormap: the colormap to use for highlighting the molecule
    :param execlude_exocyclic_from_highlight: since exocyclic atoms are not used for fc classification, this option allows the user to mask exocyclic atoms from highlights
    :param output_file: if speficied, the image is saved (format is deduced from extension)
    :param img_size: the size of the resulting Image
    :param debug: display atom indices on the structure
    :param svg: use SVG format instead of PNG
    :return: an Image of the highlighted molecule
    """
    # if no colormap is provided, do not highlight any atoms
    if colormap is None:
        colormap = ColorMap(mol, {})

    # draw
    if svg:
        d2d = rdMolDraw2D.MolDraw2DSVG(img_size[0], img_size[1])
    else:
        d2d = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    d2d.drawOptions().addAtomIndices = debug
    d2d.drawOptions().legendFontSize = 24
    d2d.drawOptions().padding = 0.1
    d2d.DrawMoleculeWithHighlights(mol, legend, colormap.atoms, colormap.bonds, {}, {})
    d2d.FinishDrawing()
    img = d2d.GetDrawingText()

    # export img
    if output_file is not None:
        output_ext = output_file.split('.')[-1].upper()
        if output_ext == 'SVG' and not svg:
            raise ValueError("Error! output file extension is SVG but image format is PNG!")
        if output_ext == 'SVG':
            with open(output_file, 'w') as SVG:
                SVG.write(img)
        else:
            raise ValueError(f"Error! Unsupported extension '{output_ext}'!")

    return img


def mols(mols: List[Mol],
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

    ..warning:: This function is based on the old RDKit drawing code (< 2020.03) and thus it is necessary to blend colors (i.e. colormap.blend()) first to represent common atoms/bonds between fragments.

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
            raise ValueError("Error! output file extension is SVG but image format is PNG!")
        if output_ext == 'SVG':
            with open(output_file, 'w') as SVG:
                SVG.write(img.data)
        else:
            raise ValueError(f"Error! Unsupported extension '{output_ext}'!")

    return img


def _get_edge_info(G: Graph, edge_attributes: List[str], attribute_names: bool, label_node_names_on_edges: bool) -> Dict:
    """
    Use the first associated data of edges for edge labelling of a networkx graph.

    :param G: a Fragment Combination graph
    :param edge_attributes: a list of edge attributes to represent on the figure
    :param attribute_names: display the attribute names on the figure (name: value)
    :return: a Dict of syntax {(node1, node2): data}
    """
    d = {}

    for edge in list(G.edges(data=True)):
        # determine what properties to keep by their names (defined as list in G attr)
        if edge_attributes is None:
            data = edge[2]
        else:
            data = {k: v for k, v in edge[2].items() if k in edge_attributes}
        # cannot use a dict for labelling edges, so just display values
        if attribute_names:
            data = [f"{k}: {v}" for k, v in data.items()]
        else:
            data = list(data.values())

        if label_node_names_on_edges:
            data.append(f"s: {edge[0]}")
            data.append(f"t: {edge[1]}")

        # format them
        d[(edge[0], edge[1])] = '; '.join(data)
    return d


def graph(G: Graph,
          colormap_nodes: List[Tuple[float]] = None,
          output_file: str = None,
          edge_attributes: List[str] = ['fcc'],
          attribute_names: bool = False,
          orientate: bool = False,
          label_node_names_on_edges: bool = False) -> Figure:
    """
    Return a matplotlib Figure of a networkx graph.

    :param G: a networkx Graph object of the fragment combinations
    :param colormap_nodes: a colormap of RGB values for the nodes (i.e. [(0, 0, 1), (0, 1, 0)]) or a ColorMap object
    :param output_file: if speficied, the image is saved (format is deduced from extension)
    :param edge_attributes: a list of edge attributes to represent on the figure
    :param attribute_names: display the attribute names on the figure (name: value)
    :return: a matplotlib Figure object
    """
    if isinstance(G, base64.bytes_types):
        G = utils.decode_object(G)

    # to orientate the graph, this is uselessly complicated, do not use this on larger networks!
    if orientate:
        H = nx.DiGraph(G, data=True)  # this creates not only node1 -> node2 but also node2 -> node1
        H.remove_edges_from(G.edges())  # this removes node1 -> node2 because it is found in the undirected graph
        G = H.reverse()  # this transforms the remaining node2 -> node1 into node1 -> node2

    if colormap_nodes is None:
        # define a 2D list instead of a single tuple to avoid matplotlib warning
        colormap_nodes = [(0.7, 0.7, 0.7)] * len(list(G.nodes()))
    elif isinstance(colormap_nodes, ColorMap):
        # define a list of colors mapped to the node iteration in G
        val_map = {k: v[0] for k, v in colormap_nodes.fragments.items()}
        colormap_nodes = [val_map.get(node, 0.0) for node in G.nodes()]  # if fragment id is not found in colormap, paint node in black instead

    pos = nx.spring_layout(G)
    edges_info = _get_edge_info(G, edge_attributes, attribute_names, label_node_names_on_edges)
    figure = plt.figure(figsize=(8, 8))
    nx.draw(G,
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
    nx.draw_networkx_edge_labels(G,
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


def rescale(mol: Mol, f: float = 1.4):
    """Rescale the coordinates of a Mol with a factor f.

    :param mol: a Mol which is modified in place.
    :param f: the factor for rescaling coordinates
    """
    tm = np.zeros((4, 4), np.double)
    for i in range(3):
        tm[i, i] = f
    tm[3, 3] = 1.0
    AllChem.TransformMol(mol, tm)


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


def reaction(mol1: Mol, mol2: Mol, sub_img_size: tuple = (200, 200), svg: bool = True, output_file: str = None):
    """Wrapper function around RDKit ReactionToImage function.
    If the molecules are Mol objects, they are converted to Smiles. If not, they are
    assumed to be already Smiles.

    .. warning:: There is currently no way of not displaying aromatic rings instead of kekulized rings. The original SMILES for the reaction can be displayed using a DEBUG logging level.

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
    logging.debug("rxn_str='%s'", rxn_str)
    rxn = rdChemReactions.ReactionFromSmarts(rxn_str, useSmiles=True)
    img = Draw.ReactionToImage(rxn, subImgSize=sub_img_size, useSVG=svg)

    # export img
    if output_file is not None:
        output_ext = output_file.split('.')[-1].upper()
        if output_ext == 'SVG' and not svg:
            raise ValueError(f"Error! output file extension is SVG but image format is PNG!")
        if output_ext == 'SVG':
            with open(output_file, 'w') as SVG:
                SVG.write(img)  # inconsistent svg attributes between mols and reactions
        else:
            raise ValueError(f"Error! Unsupported extension '{output_ext}'!")
    return img


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ColorMap:
    """A class containing all the required information for highlighting a Molecule.
    It is represented by the count of fragments, atom- and bond colors.
    """

    def __init__(self, mol: Mol, d_aidxs: dict, palette: str = None):
        """
        :param mol: the molecule to highlight. Atom/Bond properties '_color' and 'num_colors' are modified in place.
        :param d_aidxs: a dictionary containing fragment ids as keys and molecule atom indices as values, i.e. {'frag1': [(0, 1, 2)]. 'frag2': [(2, 3, 4), (5, 6, 7)]}
        :param palette: a seaborn palette defined by a string. A list of possible palette names can be found at: http://www.python-simple.com/img/img45.png. If none is provided, an intern palette is used instead. Please note that a same fragment occuring mulitple times will be displayed with only one color, but each time with a 10% darker shade. This can happen up to 5x, after what the color shade is reset to default to prevent colors from turning pitch black.
        """
        colormap_atoms = {}  # atoms
        colormap_bonds = {}  # bonds
        colormap_fragments = {}  # fragments are colored and stored in order
        if palette is None:
            palette = DEFAULT_PALETTE
        else:
            palette = sns.color_palette(palette)
            palette = {f"COL_{str(i+1).zfill(2)}": color for i, color in enumerate(palette)}
        # pick a color for each fragment type
        for i, (fragment_id, aidxs_l) in enumerate(d_aidxs.items()):
            colors_k = list(palette.keys())  # colors are sorted, so always red, then green, etc.
            color = palette[colors_k[i % len(colors_k)]]  # the great idea here is to loop back to the red once all colors have been used
            colormap_fragments[fragment_id] = []  # frag1: [(0, 0, 1)]
            # pick a shade for each occurrence of a same fragment type
            k = 0
            for j, aidxs in enumerate(aidxs_l):
                # new_color = color
                # new_color = color  # I actually find it harder to justify different shades
                new_color = tuple((x * (1.0 - 0.05 * k) for x in color))  # 5% darker for each fragment of the same type  ## TODO: define a range for j when colors get 10% lighter instead
                colormap_fragments[fragment_id].append(new_color)
                # color atoms
                for aidx in aidxs:
                    if aidx not in colormap_atoms.keys():
                        colormap_atoms[aidx] = [new_color]
                    else:
                        colormap_atoms[aidx].append(new_color)

                # color bonds
                bidxs = self._get_bidxs(mol, aidxs)
                for bidx in bidxs:
                    if bidx not in colormap_bonds.keys():
                        colormap_bonds[bidx] = [new_color]
                    else:
                        colormap_bonds[bidx].append(new_color)

                # continue to darken only if it does not end up being black
                if j < 6:
                    k += 1
                else:
                    k = 0  # reset to default

        # init attributes
        self.fragments = colormap_fragments
        self.atoms = colormap_atoms
        self.bonds = colormap_bonds
        self.palette = palette

    def __repr__(self):
        """Return a string representation of the ColorMap object
        """
        s = 'ColorMap={'
        frags = []
        # listing a huge dictionary of tuples is not that great for representing data, so just display fragment colors
        for i, v in enumerate(self.fragments.items()):
            # find the corresponding key in the palette for the color value, so one can simply display the color name
            frags = [f"{v[0]}: {list(self.palette.keys())[list(self.palette.values()).index(v[1][0])]}" for i, v in enumerate(self.fragments.items())]
        return f"{s}{', '.join(frags)}" + '}'

    def _get_bidxs(self, mol: Mol, aidxs: Set[int]) -> Set[int]:
        """From an iterable of atom indices (aidxs), return a set of bond indices (bidxs)
        found between atoms in aidxs.

        :param mol: the molecule to highlight
        :param aidxs: a set with the atom indices of the fragment to highlight
        :return: a set with the corresponding bond indices
        """
        # get all bonds attached to atom indices (as set: {0, 1, 2, 3})
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

    def blend(self):
        """Blend colors found in a ColorMap.

        The ColorMap object is modified in place.
        """

        for k, v in self.atoms.items():
            self.atoms[k] = self._blend_multiple_colors(v)
        for k, v in self.bonds.items():
            self.bonds[k] = self._blend_multiple_colors(v)

    def _blend_two_colors(self, color1: Tuple[int], color2: Tuple[int], alpha: float = 0.5) -> Tuple[int]:
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
        return tuple((sqrt((1-alpha) * (x1 ** 2) + alpha * (x2 ** 2)) for x1, x2 in zip(color1, color2)))

    def _blend_multiple_colors(self, colors: List[Tuple[int]]) -> List[Tuple[int]]:
        """Given an iterable of colors represented in RGB format (i.e. [(1, 0, 0), (0, 1, 0), (0, 0, 1)]),
        iteratively apply the blend_two_colors function.
        For instance, if colors are red, green, and blue (respectively from above), then:
            - first blend red and green applying an alpha
        """
        num_colors = len(colors)
        # if just 1 color, just return it as it is
        if num_colors < 2:
            return colors

        # iterate over colors
        color_result = colors[0]  # initiate at the first color
        for i in range(1, len(colors)):
            color_to_add = colors[i]
            old_alpha = 1 - (1/(i + 1))
            alpha = 1 - old_alpha
            color_result = self._blend_two_colors(color_result, color_to_add, alpha)
            # when I mix red and green, I want to get some yellow, not some olive color,
            # so here I hardcode the yellow I long for (golden yellow)
            if tuple(round(x, 4) for x in color_result) == (0.7211, 0.8246, 0.4472):
                color_result = (1.0, 0.9294, 0.0)
            # in case of the same color being applied, use a 10% darker shade
            if color_result == color_to_add:
                color_result = tuple((x * 0.90 for x in color_result))

        return [color_result]


class DepictionValidator:
    """
    Toolkit for estimation of depiction quality.

    This is not my code, it was copied from:
    https://gitlab.ebi.ac.uk/pdbe/ccdutils/blob/master/pdbeccdutils/core/depictions.py

    I did this because everytime I try to create an environment using the pdbeccdutils library,
    I have issues and I have to try to use an older version. So I just extracted the little piece
    that is required for my project, so I don't have problems anymore.
    """

    def __init__(self, mol):
        self.mol = mol
        self.conformer = mol.GetConformer()
        self.bonds = self.mol.GetBonds()

        atoms = [self.conformer.GetAtomPosition(i) for i in range(0, self.conformer.GetNumAtoms())]
        atom_centers = [[atom.x, atom.y, atom.z] for atom in atoms]

        self.kd_tree = KDTree(atom_centers)

    def _intersection(self, bondA, bondB):
        """
        True if two bonds collide, false otherwise. Note that false is
        retrieved even in case the bonds share common atom, as this is
        not a problem case. Cramer's rule is used for the linear
        equations system.

        Args:
            bondA (rdkit.Chem.rdchem.Bond): this bond
            bondB (rdkit.Chem.rdchem.Bond): other bond

        Returns:
            bool: true if bonds share collide, false otherwise.
        """
        atoms = [bondA.GetBeginAtom(), bondA.GetEndAtom(), bondB.GetBeginAtom(), bondB.GetEndAtom()]
        names = [a.GetProp('name') for a in atoms]
        points = [self.conformer.GetAtomPosition(a.GetIdx()) for a in atoms]

        vecA = Geometry.Point2D(points[1].x - points[0].x, points[1].y - points[0].y)
        vecB = Geometry.Point2D(points[3].x - points[2].x, points[3].y - points[2].y)

        # we need to set up directions of the vectors properly in case
        # there is a common atom. So we identify angles correctly
        # e.g. B -> A; B -> C and not A -> B; C -> B.
        if len(set(names)) == 3:
            angle = self.__get_angle(names, vecA, vecB)
            return angle < 10.0

        # Cramer's rule to identify intersection
        det = vecA.x * -vecB.y + vecA.y * vecB.x
        if round(det, 2) == 0.00:
            return False

        a = points[2].x - points[0].x
        b = points[2].y - points[0].y

        detP = (a * -vecB.y) - (b * -vecB.x)
        p = round(detP / det, 3)

        if (p < 0 or p > 1):
            return False

        detR = (vecA.x * b) - (vecA.y * a)
        r = round(detR / det, 3)

        if 0 <= r <= 1:
            return True

        return False

    def __find_element_with_max_occurrence(self, array):
        """Find element with most occurrences in the list

        Args:
            array (list of str): Array to be searched

        Returns:
            str: Value with most occurrences in the list
        """
        temp = {}
        for i in array:
            if i in temp:
                temp[i] += 1
            else:
                temp[i] = 1

        max_occur = max(temp.values())

        for k, v in temp.items():
            if v == max_occur:
                return k

        return ''

    def __get_angle(self, names, vecA, vecB):
        """Get the size of the angle formed by two bonds which share
        common atom.

        Args:
            names (list of str): List of atom names forming bonds
                [A, B, C, D] for AB and CD.
            vecA (Geometry.Point2D): Vector representing AB bond.
            vecB (Geometry.Point2D): Vector representing CD bond.

        Returns:
            float: Size of the angle in degrees.
        """
        pivot = self.__find_element_with_max_occurrence(names)

        if names[0] != pivot:  # Atoms needs to be order to pick vectors correctly
            vecA = vecA * -1

        if names[2] != pivot:
            vecB = vecB * -1

        radians = vecA.AngleTo(vecB)
        angle = 180 / math.pi * radians

        return angle

    def has_degenerated_atom_positions(self, threshold):
        """
        Detects whether the structure has a pair or atoms closer to each
        other than threshold. This can detect structures which may need
        a template as they can be handled by RDKit correctly.

        Arguments:
            threshold (float): Bottom line to use for spatial search.

        Returns:
            (bool): if such atomic pair is found
        """

        for i in range(0, len(self.conformer.GetNumAtoms())):
            center = self.conformer.GetAtomPosition(i)
            point = [center.x, center.y, center.z]
            surrounding = self.kd_tree.query_ball_point(point, threshold)

            if len(surrounding) > 1:
                return True

        return False

    def count_suboptimal_atom_positions(self, lowerBound, upperBound):
        """
        Detects whether the structure has a pair or atoms in the range
        <lowerBound, upperBound> meaning that the depiction could
        be improved.

        Arguments:
            lowerBound (float): lower bound
            upperBound (float): upper bound

        Returns:
            bool: indication whether or not the atoms are not in
            optimal coordinates
        """
        counter = 0
        for i in range(self.conformer.GetNumAtoms()):
            center = self.conformer.GetAtomPosition(i)
            point = [center.x, center.y, center.z]
            surroundingLow = self.kd_tree.query_ball_point(point, lowerBound)
            surroundingHigh = self.kd_tree.query_ball_point(point, upperBound)

            if len(surroundingHigh) - len(surroundingLow) > 0:
                counter += 1

        return counter / 2

    def count_bond_collisions(self):
        """
        Counts number of collisions among all bonds. Can be used for estimations of how 'wrong'
        the depiction is.

        Returns:
            int: number of bond collisions per molecule
        """

        errors = 0

        for i in range(0, len(self.bonds)):
            for a in range(i + 1, len(self.bonds)):
                result = self._intersection(self.bonds[i], self.bonds[a])

                if result:
                    errors += 1
        return errors

    def has_bond_crossing(self):
        """
        Tells if the structure contains collisions

        Returns:
            bool: Indication about bond collisions
        """
        return self.count_bond_collisions() > 0

    def depiction_score(self):
        """
        Calculate quality of the ligand depiction. The higher the worse.
        Ideally that should be 0.

        Returns:
            float: Penalty score.
        """

        collision_penalty = 1
        degenerated_penalty = 0.4

        bond_collisions = self.count_bond_collisions()
        degenerated_atoms = self.count_suboptimal_atom_positions(0.0, 0.5)

        score = collision_penalty * bond_collisions + degenerated_penalty * degenerated_atoms

        return round(score, 1)
