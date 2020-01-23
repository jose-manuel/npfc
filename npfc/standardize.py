"""
Module standardize
===================
This modules is used to standardize molecules and molecular DataFrames.
"""

# standard
import logging
import timeout_decorator
import copy
# data handling
import json
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdinchi
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize.metal import MetalDisconnector
from rdkit.Chem.MolStandardize.charge import Uncharger
from rdkit.Chem.MolStandardize.normalize import Normalizer
from rdkit.Chem.MolStandardize.tautomer import TautomerCanonicalizer
# dev library
from npfc import utils
from npfc.filter import Filter

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

TIMEOUT = 10  # TODO: find a way to specify a timeout value for a given Standardizer instance

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class FullUncharger(Uncharger):
    """A class derived from rdkit.Chem.MolStandardize.charge.Uncharger, so
    instead of attempting to create zwitterions all possible charges are removed
    from the molecule.

    For instance:

    >>> # Uncharger:
    >>> [O-][N+](C)(C)C[O-] -> [O-][N+](C)(C)CO

    >>> # FullUncharger:
    >>> [O-][N+](C)(C)C[O-] -> O[N+](C)(C)CO
    """

    def __init__(self):
        """Create an instance of FullUncharger.

        .. todo:: This will remove charges from -2 to +2 only. This could be improved using more general smarts?
        """
        # some smarts to use to find charges
        self.q_pos_1 = Chem.MolFromSmarts("[*;+]")
        self.q_pos_2 = Chem.MolFromSmarts("[*;+2]")
        self.q_neg_1 = Chem.MolFromSmarts("[*;-]")
        self.q_neg_2 = Chem.MolFromSmarts("[*;-2]")
        logging.debug(f"Initialized a new FullUncharger object")

    def full_uncharge(self, mol: Mol) -> Mol:
        """Neutralize molecule by adding/removing hydrogens.
        Does not attempt to preserve zwitterions.
        For now takes into account only charges of -2 and +2.

        :param mol: the input molecule
        :return: the uncharged molecule
        """
        logging.debug(f"Uncharging a molecule")
        mol = copy.deepcopy(mol)
        # Get atom ids for matches
        p = [x[0] for x in mol.GetSubstructMatches(self.q_pos_1)]   # +1
        p += [x[0] for x in mol.GetSubstructMatches(self.q_pos_2)]  # +2
        n = [x[0] for x in mol.GetSubstructMatches(self.q_neg_1)]  # -1
        n += [x[0] for x in mol.GetSubstructMatches(self.q_neg_2)]  # -2
        # remove positive charges
        for atom in [mol.GetAtomWithIdx(x) for x in p]:
            # Remove hydrogen and reduce formal change until neutral or no more hydrogens
            while atom.GetFormalCharge() > 0 and atom.GetNumExplicitHs() > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        # remove negative charges
        for atom in [mol.GetAtomWithIdx(x) for x in n]:
            # Add hydrogen and increase formal change until neutral
            while atom.GetFormalCharge() < 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        return mol


class Standardizer(Filter):
    """A class for standardizing molecular structures. The standardization itself is based
    on a protocol that the user can modify.

    By default this protocol consists in 12 tasks applied to each molecule invidually:

        1) **initiate_mol**: check if the molecule passed the RDKit conversion
        2) **filter_empty**: filter molecules with empty structures
        3) **disconnect_metal**: break bonds involving metallic atoms, resulting in potentially several molecules per structure.
        4) **keep_best**: retrieve only the "best" molecule from a mixture, which might not always be the largest one.
        5) **filter_hac**: filter molecules with a heavy atom count not in the accepted range. By default: hac > 3.
        6) **filter_molweight**: filter molecules with a molecular weight not in the accepted range. By default: molweight <= 1000.0.
        7) **filter_nrings**: filter molecules with a number of rings (Smallest Sets of Smallest Rings or SSSR) not in the accepted range. By default: nrings > 0.
        8) **filter_medchem**: filter molecules with elements not considered as medchem. By default: elements in H, B, C, N, O, F, P, S, Cl, Br, I.
        9) **remove_isotopes**: set all atoms to their most common isotope (i.e. 14C becomes 12C which is C).
        10) **normalize**: always write the same functional groups in the same manner.
        11) **uncharge**: remove all charges on a molecule when it is possible. This is different from rdkit.Chem.MolStandardize.charge module as there is no attempt for reaching the zwitterion.
        12) **canonicalize**: enumerate the canonical tautomer.
        13) **remove_stereo**: remove all remaining stereochemistry flags on the molecule.

    Finally, duplicate entries can be filtered using InChiKeys when postprocessing the DataFrame with molecules (argument in the run_df method).

    This results in new columns in the input DataFrame:

        - the 'mol' column: updated structure (only for the protocol)
        - the 'status' column: either passed, filtered or error.
        - the 'task' column: the latest task that was applied to the molecule.

    Optionally, it is possible to compute new 2D coordinates for depicting molecules as final step of the standardization.
    To achieve this, a scoring function is assigned to each depiction to retrieve the "best" one (least amount of overlapping atoms/bonds):

        - Input (if any)
        - rdDepictor
        - Avalon
        - rdCoordgen

    A column "depiction_method" is appended to the output DataFrame so the user knows which method was considered the best.

    .. note:: For now, the user can only change the task order and the values of filters, but it would be relatively easy to add more functionality.

    .. note:: For now structures are deglycosylated before being standardized using a KIME workflow. This is far from being ideal, at the node itself contains bugs and its execution is relatively slow.
    """

    def __init__(self,
                 protocol: str = None,
                 col_mol: str = 'mol',
                 col_id: str = 'idm',
                 filter_duplicates: bool = True,
                 compute_2D: bool = False,
                 ref_file: str = None,
                 on: str = 'inchikey',
                 # ref_dataset=None,
                 suffix: str = None,
                 elements_medchem: set = {'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'}
                 ):
        """Create a Standardizer object."""
        # filter
        super(Standardizer, self).__init__()
        # standardizer
        self._elements_medchem = elements_medchem
        self._col_id = col_id
        self._col_mol = col_mol
        self._filter_duplicates = filter_duplicates
        self._compute_2D = compute_2D
        self._on = on
        self._suffix = suffix
        self._ref_file = ref_file
        self._default_protocol = {'tasks': ['filter_empty',
                                            'disconnect_metal',
                                            'keep_best',
                                            'filter_hac',
                                            'filter_molweight',
                                            'filter_nrings',
                                            'filter_medchem',
                                            'remove_isotopes',
                                            'normalize',
                                            'uncharge',
                                            'canonicalize',
                                            'remove_stereo',
                                            ],
                                  'filter_hac': 'hac > 3',
                                  'filter_molweight': 'molweight <= 1000.0',
                                  'filter_nrings': 'nrings > 0',
                                  'filter_medchem': f'elements in {", ".join(str(x) for x in self.elements_medchem)}',
                                  }
        if protocol is None:
            self._protocol = self._default_protocol
        else:
            self._protocol = protocol
        # workers
        self.metal_disconnector = MetalDisconnector()
        self.normalizer = Normalizer()
        self.full_uncharger = FullUncharger()
        self.canonicalizer = TautomerCanonicalizer()

    @property
    def protocol(self):
        return self._protocol

    @protocol.setter
    def protocol(self, protocol: str):
        # input is a json file => convert it to a dict
        if isinstance(protocol, str):
            utils.check_arg_config_file(protocol)
            with open(protocol) as f:
                protocol = json.load(f)
        # input is a dict
        if isinstance(protocol, dict):
            if 'tasks' not in protocol.keys():
                raise ValueError("invalid protocol format (no 'tasks' key found)")
            elif not isinstance(protocol['tasks'], list) and not isinstance(protocol['tasks'], tuple):
                raise ValueError("invalid protocol format ('tasks' key is neither list or tuple)")
        # update default protocol
        self._protocol.update(protocol)

    @property
    def col_id(self) -> str:
        return self._col_id

    @col_id.setter
    def col_id(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! col_id cannot be '{value}'.")
        self._col_id = value

    @property
    def col_mol(self) -> str:
        return self._col_mol

    @col_mol.setter
    def col_mol(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! col_mol cannot be '{value}'.")
        self._col_mol = value

    @property
    def elements_medchem(self) -> set:
        return self._elements_medchem

    @elements_medchem.setter
    def elements_medchem(self, value: set) -> None:
        if not isinstance(value, set):
            raise ValueError(f"Error! elements_medchem should be a set of strings, not '{value}' ({type(value)}).")
        elif not all([isinstance(v, str) for v in value]):
            raise ValueError(f"Error! elements_medchem should be a set of strings, not '{value}' ({type(value)}).")
        self._elements_medchem = value

    @property
    def filter_duplicates(self) -> bool:
        return self._filter_duplicates

    @filter_duplicates.setter
    def filter_duplicates(self, value: bool) -> None:
        utils.check_arg_bool(value)
        self._filter_duplicates = value

    @property
    def compute_2D(self) -> bool:
        return self._compute_2D

    @compute_2D.setter
    def compute_2D(self, value: bool) -> None:
        utils.check_arg_bool(value)
        self._compute_2D = value

    @property
    def on(self) -> str:
        return self._on

    @on.setter
    def on(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! on cannot be '{value}'.")
        self._on = value

    @property
    def suffix(self) -> str:
        return self._suffix

    @suffix.setter
    def suffix(self, value: str) -> None:
        if value is not None and not isinstance(value, str):
            raise ValueError(f"Error! Either None or a str are expected for suffix, not '{value}' ({type(value)}).")
        self._suffix = value

    @property
    def ref_file(self) -> str:
        return self._ref_file

    @ref_file.setter
    def ref_file(self, value: str) -> None:
        if value is not None and not isinstance(value, str):
            raise ValueError(f"Error! Either None or a str are expected for ref_file, not '{value}' ({type(value)}).")
        self._ref_file = value

    def remove_isotopes(self, mol: Mol) -> Mol:
        """Return a molecule without any isotopes.

        :param mol: the input molecule
        :return: the molecule without isotope
        """
        mol = Mol(mol)
        for a in mol.GetAtoms():
            a.SetIsotope(0)
        return mol

    def keep_best(self, mol: Mol) -> Mol:
        """Return the "best" molecule found in a molecular structure.

        The "best" molecule is determined by the following criteria, sorted by priority:

            1) contains only medchem elements
            2) contains at least one ring
            3) has the largest molecular weight of the mixture

        To summarize:

        .. math::
            medchem > non linear > molecular weight

        So the largest molecule of a mixture might not always be selected, for instance
        a very long aliphatic chain would be dismissed to keep a benzene molecule instead.

        This is implemented in such a way because our fragments used for substructure search contain at least one ring.
        On the contrary, this long aliphatic chain would be kept in a mixture with a non-medchem molecule.

        :param mol: the input molecule(s)
        :return: the best molecule
        """
        submols = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        # no need to look further if we have only one submol!
        if len(submols) < 2:
            return mol
        # otherwise, we have to compare the submols
        # init
        logging.debug(f"found {len(submols)} submols")
        best_molweight = -1.0  # so we are sure to update this on the first iteration
        best_submol = None
        best_is_medchem = False
        best_is_non_linear = False
        # begin
        for i, submol in enumerate(submols):
            # is_medchem
            is_medchem = self.filter_mol(submol, f'elements in {", ".join(str(x) for x in self.elements_medchem)}')
            is_non_linear = self.filter_mol(submol, f"nrings > 0")
            # molweight
            molweight = Descriptors.ExactMolWt(submol)
            logging.debug(f"submol #{i}: IM={is_medchem} INL={is_non_linear} MW: {molweight}")
            # compare to the current best fragment
            update_best = False
            compute_diff = False  # check which
            # 2 criteria more important than molecular weight: is_medchem > is_non_linear
            if not best_is_medchem and is_medchem:
                update_best = True
            elif best_is_medchem and not is_medchem:
                continue
            elif not best_is_medchem and not is_medchem:
                if not best_is_non_linear and is_non_linear:
                    update_best = True
                elif best_is_non_linear and not is_non_linear:
                    continue
                else:
                    compute_diff = True
            else:  # best_is_medchem and is_medchem
                if not best_is_non_linear and is_non_linear:
                    update_best = True
                elif best_is_non_linear and not is_non_linear:
                    continue
                else:
                    compute_diff = True

            # check molweights only in case of doubt
            if not update_best and compute_diff and molweight > best_molweight:
                update_best = True

            # update best with the properties of the current mol
            if update_best:
                best_is_medchem = is_medchem
                best_is_non_linear = is_non_linear
                best_submol = submol
                best_molweight = molweight

        return best_submol

    @timeout_decorator.timeout(TIMEOUT)
    def _run(self, mol: Mol) -> tuple:
        """Helper function for run.
        Contains all tasks defined within the protocol. Since some operations are
        expansive and could last a very long time for complex molecules (normalize, canonicalize),
        a timeout value is set globally. The run function is the one that can catch the exception
        raised by timeouts.

        :param mol: the input molecule
        :return: a tuple containing the molecule, its status and the further task name it reached
        """

        # initiate_mol
        if mol is None:
            return (mol, 'error', 'initiate_mol')

        # begin protocol
        mol = copy.deepcopy(mol)
        for task in self._protocol['tasks']:
            # filter_empty
            if task == 'filter_empty':
                try:
                    if not mol.GetNumAtoms():
                        return ('', 'filtered', task)
                except ValueError:
                    return ('', 'error', task)

            # disconnect_metal
            elif task == 'disconnect_metal':
                try:
                    mol = self.metal_disconnector.disconnect(mol)
                except ValueError:
                    return (mol, 'error', 'disconnect_metal')

            # keep_best
            elif task == 'keep_best':
                try:
                    mol = self.keep_best(mol)
                except ValueError:
                    return (mol, 'error', task)

            # filters
            elif task == 'filter_medchem' or task == 'filter_hac' or task == 'filter_molweight' or task == 'filter_nrings':
                try:
                    if not self.filter_mol(mol, self._protocol[task]):
                        return (mol, 'filtered', task)
                except ValueError:
                    return (mol, 'error', task)

            # sanitize
            elif task == 'sanitize':
                try:
                    Chem.SanitizeMol(mol)
                except ValueError:
                    return (mol, 'error', task)

            # remove_isotopes
            elif task == 'remove_isotopes':
                try:
                    mol = self.remove_isotopes(mol)
                except ValueError:
                    return (mol, 'error', task)

            # normalize
            elif task == 'normalize':
                try:
                    mol = self.normalizer.normalize(mol)
                except ValueError:
                    return (mol, 'error', task)

            # uncharge
            elif task == 'uncharge':
                try:
                    mol = self.full_uncharger.full_uncharge(mol)
                except ValueError:
                    return (mol, 'error', task)

            # canonicalize
            elif task == 'canonicalize':
                # canonicalize
                try:
                    mol = self.canonicalizer.canonicalize(mol)
                except (ValueError, RuntimeError):
                    return (mol, 'error', task)

            # remove_stereo
            elif task == 'remove_stereo':
                try:
                    rdmolops.RemoveStereochemistry(mol)
                except ValueError:
                    return (mol, 'error', task)

            # something else?
            else:
                raise ValueError(f"Unknown task: {task}")

        # a molecule that passed all the protocole!
        return (mol, 'passed', 'standardize')

    def run(self, mol: Mol) -> tuple:
        """Execute the standardization protocol on a molecule.
        Molecule that exceed the timeout value are filtered with a task='timeout'.

        As a final step of the protocol, InChiKeys ('inchikey') are computed for identifying molecules.

        :param mol: the input molecule
        :return: a tuple containing the molecule, its status and the further task name it reached
        """
        try:
            return self._run(mol)
        except timeout_decorator.TimeoutError:
            return (mol, 'filtered', 'timeout')

    def run_df(self, df: DataFrame) -> tuple:
        """Apply the standardization protocol on a DataFrame, with the possibility of directly filtering duplicate entries as well.
        This can be very useful as the standardization process can expose duplicate entries due to salts removal, neutralization,
        canonical tautomer enumeration, and stereochemistry centers unlabelling

        If a reference file is specified, duplicate removals becomes possible accross chunks.


        :param df: the input DataFrame
        :return: three DataFrames separated by status:

            - passed
            - filtered
            - error

        .. note:: As a side effect, the output DataFrames get indexed by idm. The 'inchikey' col is not returned, but the values can be accessed using the reference file.
        """
        # run standardization protocol
        df.index = df[self.col_id]
        df.loc[:, self.col_mol], df.loc[:, 'status'], df.loc[:, 'task'] = zip(*df[self.col_mol].map(self.run))
        # do not apply filter duplicates on molecules with errors or that were already filtered for x reasons
        df_error = df[df['status'] == 'error']
        df_filtered = df[df['status'] == 'filtered']
        df = df[df['status'] == 'passed']
        df = df.copy()  # only way I found to suppress pandas warnings in a "clean" way

        # compute InChiKeys
        df.loc[:, 'inchikey'] = df.loc[:, self.col_mol].map(rdinchi.MolToInchiKey)

        # tuple of dataframes
        return (df, df_filtered, df_error)
