"""
Module test_06_draw
====================
Tests for drawing molecules.
"""
# standard
from pathlib import Path
import warnings
import filecmp
# data handling
from collections import OrderedDict
import pickle
import base64
import pandas as pd
# chemoinformatics
from rdkit import Chem
from rdkit import RDLogger
# tests
import pytest
from npfc.load import decode_mol_base64
from npfc import draw
# configure logging
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
# debug
# import logging
# logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def ref_file_highlight_mol():
    return "npfc/data/test_06_draw_highlight_mol.svg"


@pytest.fixture
def df_map():
    df_map = pd.read_csv("tests/test_case_chembl_2_map.csv.gz", compression="gzip", sep="|")
    df_map["mol"] = df_map["mol"].map(decode_mol_base64)
    df_map["graph"] = df_map["graph"].map(lambda x: pickle.loads(base64.b64decode(x)))
    df_map["colormap"] = df_map["colormap"].map(lambda x: pickle.loads(base64.b64decode(x)))

    return df_map


@pytest.fixture
def mols_2D():
    return {
            "simple": Chem.MolFromSmiles("C1CCCCCNc2cc[n+](Cc3ccc(Cc4ccc(C[n+]5ccc(NCCCC1)c6ccccc56)cc4)cc3)c7ccccc27"),  # CHEMBL446445
            "medium": Chem.MolFromSmiles("CC1(C)O[C@@H]2[C@]3(CCC[C@@H]4[C@H]5CO[C@H]6C=C(CO)[C@@]2(O)[C@@]34[C@@H]56)O1"),  # CHEMBL1269998
            "hard": Chem.MolFromSmiles("C[C@]12CC[C@@H]3OC(=O)[C@@]45C[C@@H](C[C@@H](O)[C@H]4[C@]36CO[C@H](OC1)[C@H]26)C7=C5O[C@@]8(CC7)[C@@H]9C[C@@H](O)[C@H]%10[C@@]%11%12CO[C@@H](O)[C@@H]%11C(C)(C)[C@@H](O)C[C@@H]%12OC(=O)[C@]%10(C9)C8=O"),  # CHEMBL1079125
            }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_compute_2D(mols_2D):
    """Test if given methods yield the best 2D representations of molecules. Random seeds are not fixed but a run of 100x was applied with no variation in results."""
    # simple case where CoordGen does the trick just fine
    mol_s = draw.compute_2D(mols_2D["simple"])
    assert mol_s.GetProp("_2D") == "CoordGen"
    # medium case where CoordGen actually performs worse than Avalon
    mol_m = draw.compute_2D(mols_2D["medium"])
    assert mol_m.GetProp("_2D") == "Avalon"
    # complex case where none of the tools yields a perfect score, but CoordGen is the "least worse"
    mol_c = draw.compute_2D(mols_2D["hard"])
    assert mol_c.GetProp("_2D") == "rdDepictor"


def test_highlight_mol(df_map, ref_file_highlight_mol):
    """Highlight a molecule using the colormap object defined by the previous test file. Compare the output SVG to a reference SVG."""
    test_file_highlight_mol = "tests/tmp/test_06_draw_highlight_mol.svg"
    # check if ref file is available
    assert Path(ref_file_highlight_mol).exists() is True
    # check colormap
    row = df_map.iloc[0]
    colormap = row["colormap"]
    assert colormap.fragments == OrderedDict([('1141', (1.0, 0.6, 0.6)), ('678', (0.2, 1.0, 0.2))])
    assert colormap.atoms == {2: (1.0, 0.6, 0.6),
                              3: (1.0, 0.6, 0.6),
                              4: (1.0, 0.6, 0.6),
                              5: (1.0, 0.6, 0.6),
                              6: (1.0, 0.6, 0.6),
                              7: (1.0, 0.6, 0.6),
                              8: (1.0, 0.9294, 0.0),
                              9: (0.2, 1.0, 0.2),
                              10: (0.2, 1.0, 0.2),
                              21: (0.2, 1.0, 0.2),
                              22: (0.2, 1.0, 0.2),
                              23: (0.2, 1.0, 0.2),
                              24: (1.0, 0.9294, 0.0),
                              25: (1.0, 0.6, 0.6),
                              26: (1.0, 0.6, 0.6),
                              27: (1.0, 0.6, 0.6),
                              28: (1.0, 0.9294, 0.0),
                              29: (1.0, 0.6, 0.6),
                              30: (1.0, 0.6, 0.6),
                              31: (1.0, 0.6, 0.6),
                              }
    assert colormap.bonds == {2: (1.0, 0.6, 0.6),
                              3: (1.0, 0.6, 0.6),
                              4: (1.0, 0.6, 0.6),
                              5: (1.0, 0.6, 0.6),
                              6: (1.0, 0.6, 0.6),
                              7: (1.0, 0.6, 0.6),
                              8: (0.2, 1.0, 0.2),
                              9: (0.2, 1.0, 0.2),
                              20: (0.2, 1.0, 0.2),
                              21: (0.2, 1.0, 0.2),
                              22: (0.2, 1.0, 0.2),
                              23: (0.2, 1.0, 0.2),
                              24: (1.0, 0.6, 0.6),
                              25: (1.0, 0.6, 0.6),
                              26: (1.0, 0.6, 0.6),
                              27: (1.0, 0.6, 0.6),
                              28: (1.0, 0.6, 0.6),
                              29: (1.0, 0.6, 0.6),
                              30: (1.0, 0.6, 0.6),
                              32: (1.0, 0.6, 0.6),
                              33: (1.0, 0.6, 0.6),
                              34: (1.0, 0.6, 0.6),
                              35: (1.0, 0.9294, 0.0),
                              36: (0.2, 1.0, 0.2),
                              37: (1.0, 0.9294, 0.0),
                              0: (1, 1, 1),
                              1: (1, 1, 1),
                              38: (1, 1, 1),
                              10: (1, 1, 1),
                              11: (1, 1, 1),
                              12: (1, 1, 1),
                              13: (1, 1, 1),
                              14: (1, 1, 1),
                              15: (1, 1, 1),
                              16: (1, 1, 1),
                              17: (1, 1, 1),
                              18: (1, 1, 1),
                              19: (1, 1, 1),
                              31: (1, 1, 1),
                              }
    # compute a SVG drawing of the molecule with highlights
    img = draw.highlight_mol(row["mol"], colormap=colormap, svg=True)
    with open(test_file_highlight_mol, "w") as FILE:
        FILE.write(img)
    # check if generated file is identical to ref file (size and content)
    assert filecmp.cmp(test_file_highlight_mol, ref_file_highlight_mol) is True
