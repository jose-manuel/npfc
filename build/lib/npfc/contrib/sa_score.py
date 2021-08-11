"""
Module sa_score
=================
This modules contains the functions provided by Ertl et al at:

https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py

It was simply edited to suit my needs and tastes, I do not take credit for it.

It uses a pre-computed model to score structures based on the presence or
absence of fragments.
"""
# standard
import math
import gzip
# data handling
import pickle
# chemoinformatics
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdMolDescriptors


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class SAScorer:
    """A class for computing the Synthetic Accessibility score, after Ertl et al.
    Neither the algorithm nor the model are mine, I just downloaded them from the ref
    and adapted the code for my project.

    The original model is file is required for scoring molecules, its location
    needs to be provided as argument.

    Reference:
    https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py
    """

    def __init__(self, file_model: str):
        self.fscores = self.load_model(file_model)
        self.file_model = file_model

    def __repr__(self):
        return f"SAScorer with model loaded at '{self.file_model}'"

    def load_model(self, file_model: str):
        data = pickle.load(gzip.open(file_model))
        outDict = {}
        for i in data:
            for j in range(1, len(i)):
                outDict[i[j]] = float(i[0])
        return outDict

    def _numBridgeheadsAndSpiro(self, mol, ri=None):
        nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        return nBridgehead, nSpiro

    def score(self, m: Mol) -> float:
        """Compute the Synthetic Accessiblity score of a molecule.
        It ranges from 1 (very easy to make) to 10 (very hard to make).

        :param mol: the input molecule
        :return: the score
        """
        # fragment score
        fp = rdMolDescriptors.GetMorganFingerprint(m, 2)  # circular fingerprint *radius*
        fps = fp.GetNonzeroElements()
        score1 = 0.
        nf = 0
        for bitId, v in fps.items():
            nf += v
            sfp = bitId
            score1 += self.fscores.get(sfp, -4) * v
        score1 /= nf

        # features score
        nAtoms = m.GetNumAtoms()
        nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
        ri = m.GetRingInfo()
        nBridgeheads, nSpiro = self._numBridgeheadsAndSpiro(m, ri)
        nMacrocycles = 0
        for x in ri.AtomRings():
            if len(x) > 8:
                nMacrocycles += 1

        sizePenalty = nAtoms**1.005 - nAtoms
        stereoPenalty = math.log10(nChiralCenters + 1)
        spiroPenalty = math.log10(nSpiro + 1)
        bridgePenalty = math.log10(nBridgeheads + 1)
        macrocyclePenalty = 0.
        # ---------------------------------------
        # This differs from the paper, which defines:
        #  macrocyclePenalty = math.log10(nMacrocycles+1)
        # This form generates better results when 2 or more macrocycles are present
        if nMacrocycles > 0:
            macrocyclePenalty = math.log10(2)

        score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

        # correction for the fingerprint density
        # not in the original publication, added in version 1.1
        # to make highly symmetrical molecules easier to synthetise
        score3 = 0.
        if nAtoms > len(fps):
            score3 = math.log(float(nAtoms) / len(fps)) * .5

        sascore = score1 + score2 + score3

        # need to transform "raw" value into scale between 1 and 10
        min = -4.0
        max = 2.5
        sascore = 11. - (sascore - min + 1) / (max - min) * 9.
        # smooth the 10-end
        if sascore > 8.:
            sascore = 8. + math.log(sascore + 1. - 9.)
        if sascore > 10.:
            sascore = 10.0
        elif sascore < 1.:
            sascore = 1.0

        return sascore
