"""
Module save
================
"""


class Saver:

    def __init__(self,
                 shuffle: bool = False,
                 random_seed: int = None,
                 chunk_size: int = None,
                 encode_mols: bool = True,
                 col_mol: str = 'mol'):
        self._shuffle = shuffle
        self._random_seed = random_seed
        self._chunk_size = chunk_size
        self._encode_mols = encode_mols
        self._col_mol = col_mol
