"""
Module np_score
=================
This modules contains the functions provided by Ertl et al at:

https://github.com/rdkit/rdkit/blob/master/Contrib/NP_Score/npscorer.py

It was simply edited to suit my needs and tastes, I do not take credit for it.

It uses a pre-computed model to score structures based on the presence or
absence of fragments. Scores vary between -5 (synthetic) and +5 (natural product).
"""
# standard
import math
import gzip
# data handling
import pickle
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import rdMolDescriptors
# docs
from typing import Union

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class NPScorer:
    """A class for computing the Natural-Product-Likeness score, after Ertl et al.
    Neither the algorithm nor the model are mine, I just downloaded them from the ref
    and adapted the code for my project.

    The original model is file is required for scoring molecules, its location
    needs to be provided as argument.

    Reference:
    https://github.com/rdkit/rdkit/blob/master/Contrib/NP_Score/npscorer.py
    """

    def __init__(self, file_model: str):
        self.fscores = pickle.load(gzip.open(file_model))
        self.file_model = file_model

    def __repr__(self):
        return f"NPScorer with model loaded at '{self.file_model}'"

    def score(self, mol: Mol, confidence=False) -> Union[float, dict]:
        """Compute the Natural-Product-Likeness of a molecule.
        This scores varies between -5.0 (synthetic) and +5.0 (natural). It relies
        on the predefined fragment scores from the model, therefore a confidence
        value (0-1.0) can be computed to estimate the proportion of the molecule that
        was actually assessed and thus, how reliable the score is.

        :param mol: the input molecule
        :param confidence: compute the confidence value
        :return: either directly the 'np_score' or a dictionary with the 'np_score' and the 'confidence' values.
        """

        if mol is None:
            raise ValueError('Error! Input molecule is None!')

        # compute fingerprint bits
        fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
        bits = fp.GetNonzeroElements()

        # calculating the score
        score = 0.0
        bits_found = 0
        for bit in bits:
            if bit in self.fscores:
                bits_found += 1
                score += self.fscores[bit]
        score /= float(mol.GetNumAtoms())

        # preventing score explosion for exotic molecules
        if score > 4:
            score = 4. + math.log10(score - 4. + 1.)
        elif score < -4:
            score = -4. - math.log10(-4. - score + 1.)

        # avoiding unnecessary float precision
        score = round(score, 4)

        if confidence:
            # confidence score to asses how much of the tested fragments of the molecule
            # represent of the molecule (the higher, the better)
            confidence = float(bits_found / len(bits))

            return {'np_score': score, 'confidence': round(confidence, 4)}
        else:
            return score
