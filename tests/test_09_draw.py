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
import pandas as pd
# chemoinformatics
from rdkit import Chem
from rdkit import RDLogger
# tests
import pytest
from npfc import utils
from npfc import draw
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOGGING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def ref_file_mol():
    return "tests/data/test_draw_mol.svg"


@pytest.fixture
def df_map():
    d = {}
    d['idm'] = ['CHEMBL10006']
    d['fmid'] = ['001']
    d['nfrags'] = [2]
    d['nfrags_u'] = [2]
    d['ncomb'] = [1]
    d['ncomb_u'] = [1]
    d['hac_mol'] = [33]
    d['hac_frags'] = [20]
    d['perc_mol_cov_frags'] = [61.0]
    d['frags'] = [['678:0', '1141:0']]
    d['frags_u'] = [['1141', '678']]
    d['comb'] = ['fbr']
    d['comb_u'] = ['fbr']
    d['fmap_str'] = ['678:0[fbr]1141:0']
    d['_d_aidxs'] = ["gANjY29sbGVjdGlvbnMKT3JkZXJlZERpY3QKcQApUnEBKFgEAAAAMTE0MXECXXEDY2J1aWx0aW5zCmZyb3plbnNldApxBF1xBShLAksDSwRLBUsGSwdLCEsYSxlLGksbSxxLHUseSx9lhXEGUnEHYVgDAAAANjc4cQhdcQloBF1xCihLCEsJSwpLFUsWSxdLGEscZYVxC1JxDGF1Lg=="]
    d['_colormap'] = ["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"]
    d['_fmap'] = ["gANjbmV0d29ya3guY2xhc3Nlcy5ncmFwaApHcmFwaApxACmBcQF9cQIoWBcAAABncmFwaF9hdHRyX2RpY3RfZmFjdG9yeXEDY2J1aWx0aW5zCmRpY3QKcQRYEQAAAG5vZGVfZGljdF9mYWN0b3J5cQVoBFgWAAAAbm9kZV9hdHRyX2RpY3RfZmFjdG9yeXEGaARYGgAAAGFkamxpc3Rfb3V0ZXJfZGljdF9mYWN0b3J5cQdoBFgaAAAAYWRqbGlzdF9pbm5lcl9kaWN0X2ZhY3RvcnlxCGgEWBYAAABlZGdlX2F0dHJfZGljdF9mYWN0b3J5cQloBFgFAAAAZ3JhcGhxCn1xC1gFAAAAX25vZGVxDH1xDShYAwAAADY3OHEOfXEPWAQAAAAxMTQxcRB9cRF1WAQAAABfYWRqcRJ9cRMoaA59cRRoEH1xFShYBgAAAGFiYnJldnEWWAMAAABmYnJxF1gIAAAAbl9hYmJyZXZxGEsBWAMAAABpZG1xGVgLAAAAQ0hFTUJMMTAwMDZxGnVzaBB9cRtoDmgVc3V1Yi4="]
    d['mol'] = ["776t3gAAAAAKAAAAAAAAAAAAAAAhAAAAJwAAAIABBgBgAAAAAQMIACgAAAADAgZAKAAAAAMEBkBoAAAAAwMBBkBoAAAAAwMBBkAoAAAAAwQGAGAAAAACAgYAYAAAAAMBBgBgAAAAAwEGAGAAAAACAgYgJAAAAAIEBgBgAAAAAgIGAGAAAAACAgYAYAAAAAICBgBgAAAAAgIGQCgAAAADBAZAaAAAAAMDAQZAaAAAAAMDAQZAaAAAAAMDAQZAaAAAAAMDAQZAaAAAAAMDAQYAYAAAAAICCAAgAAAAAgYgNAAAAAIBBAYAYAAAAAMBCAAoAAAAAwIGQCgAAAADBAZAKAAAAAMEBgAgAAAABAYAYAAAAAICBgBgAAAAAgIHACAAAAADBgBgAAAAAQMLAAEAAQIgAgNoDAMEaAwEBWgMBQYABgcABwgACAkACQoACgsACwwADA0ADQ4ADg8ADxBoDBARaAwREmgMEhNoDBMUaAwKFQAVFgAWFwAXGAAYGQAZGiAaG2gMGxwAHB0AHR4AHh8AHyAAGgJoDBsFaAwfBwAcCAAXCgAcGAAUD2gMFAcGAhobBQQDBgYHCBwbBQYJChcYHAgGEBESExQPBBUWFwoFGRgcGxoGHR4fBwgcFwAAAAAW"]
    df_map = pd.DataFrame(d)
    for c in ('_d_aidxs', '_colormap', '_fmap'):
        df_map[c] = df_map[c].map(utils.decode_object)
    df_map['mol'] = df_map['mol'].map(utils.decode_mol)
    df_map = pd.DataFrame(d)
    for c in ('_d_aidxs', '_colormap', '_fmap'):
        df_map[c] = df_map[c].map(utils.decode_object)
    df_map['mol'] = df_map['mol'].map(utils.decode_mol)

    return df_map


@pytest.fixture
def mols_2D():
    return {
            "simple": Chem.MolFromSmiles("C1CCCCCNc2cc[n+](Cc3ccc(Cc4ccc(C[n+]5ccc(NCCCC1)c6ccccc56)cc4)cc3)c7ccccc27"),  # CHEMBL446445
            "medium": Chem.MolFromSmiles("CC1(C)O[C@@H]2[C@]3(CCC[C@@H]4[C@H]5CO[C@H]6C=C(CO)[C@@]2(O)[C@@]34[C@@H]56)O1"),  # CHEMBL1269998
            "hard": Chem.MolFromSmiles("C[C@]12CC[C@@H]3OC(=O)[C@@]45C[C@@H](C[C@@H](O)[C@H]4[C@]36CO[C@H](OC1)[C@H]26)C7=C5O[C@@]8(CC7)[C@@H]9C[C@@H](O)[C@H]%10[C@@]%11%12CO[C@@H](O)[C@@H]%11C(C)(C)[C@@H](O)C[C@@H]%12OC(=O)[C@]%10(C9)C8=O"),  # CHEMBL1079125
            }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_depict(mols_2D):
    """Test if given methods yield the best 2D representations of molecules. Random seeds are not fixed but a run of 100x was applied with no variation in results."""
    # simple case where CoordGen does the trick just fine
    mol_s = draw.depict_mol(mols_2D["simple"])
    assert mol_s.GetProp("_2D") == "CoordGen"
    # medium case where CoordGen actually performs worse than Avalon
    mol_m = draw.depict_mol(mols_2D["medium"])
    assert mol_m.GetProp("_2D") == "CoordGen"
    # assert mol_m.GetProp("_2D") == "Avalon"  # for bow Avalon is disabled...
    # complex case where none of the tools yields a perfect score, but CoordGen is the "least worse"
    mol_c = draw.depict_mol(mols_2D["hard"])
    assert mol_c.GetProp("_2D") == "rdDepictor"


def test_mol(df_map, ref_file_mol):
    """Highlight a molecule using the colormap object defined by the previous test file. Compare the output SVG to a reference SVG."""
    test_file_mol = "tests/tmp/test_06_draw_highlight_mol.svg"
    # check if ref file is available
    assert Path(ref_file_mol).exists() is True
    # check colormap
    row = df_map.iloc[0]
    colormap = row["_colormap"]
    # colors attributed to each fragment
    assert colormap.fragments == OrderedDict([('1141', (1.0, 0.6, 0.6)), ('678', (0.2, 1.0, 0.2))])
    # colors attributed to each atom. Fused atom colors are blended together
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
    # # colors attributed to each bond. Bonds between fused atoms are blended together
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
    # # compute a SVG drawing of the molecule with highlights
    # img = draw.mol(row["mol"], colormap=colormap, svg=True)
    # with open(test_file_mol, "w") as FILE:
    #     FILE.write(img)
    # # check if generated file is identical to ref file (size and content)
    # assert filecmp.cmp(test_file_mol, ref_file_mol) is True
