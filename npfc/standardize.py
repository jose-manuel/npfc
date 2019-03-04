"""
Module standardize
===================
This modules is used to standardize molecules and molecular DataFrames.
"""

# standard
import logging
import timeout_decorator
import copy
from pathlib import Path
# data handling
import json
from itertools import chain
import pandas as pd
# chemoinformatics
from rdkit import Chem
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

    def full_uncharge(self, mol):
        """Neutralize molecule by adding/removing hydrogens.
        Does not attempt to preserve zwitterions.
        For now takes into account only charges of -2 and +2.

        :param mol: The input molecule.
        :type mol: : rdkit.Chem.Mol
        :return mol: The uncharged molecule.
        :rtype mol: : rdkit.Chem.Mol
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


class DuplicateFilter:
    """A class used for filtering duplicate molecules in on or several DataFrame(s)."""

    def __init__(self, on='inchikey', ref_file=None, col_mol='mol', col_id='idm'):
        """Create an instance of DuplicateFilter with following parameters:

        :param on: The property to use for identifying duplicate molecules.
        :type on: str
        :param ref_file: The path to the reference file used to store synonyms.
        :type ref_file: str
        :param col_mol: The DataFrame column name where molecules are stored.
        :type col_mol: str
        :param col_id: The DataFrame column name where molecule identifiers are stored.
        :type col_id: str

        .. note:: For now, only the 'inchikey' property works for the 'on' parameter. Moreover InchIKey are computed wether provided or not.

        """
        self._on = on
        self.col_mol = col_mol
        self.col_id = col_id
        self.col_id_synonyms = self.col_id + "_synonyms"
        self.ref_file = ref_file
        logging.debug(f"Initialized a new DuplicateFilter object")

    @property
    def on(self):
        return self._on

    @on.setter
    def on(self, value):
        if value != 'inchikey':
            raise ValueError(f"for now argument 'on' only accepts 'inchikey' as value, not '{value}'")
        self._on = value

    def find_synonyms(self, df):
        """Find synonyms within a DataFrame.

        The Synonyms DataFrame follows the follow syntax:
        - rowid: The InChIKey of the molecule ('inchikey').
        - idm_synonyms: The list of all molecule ids sharing the same InChIKey.
        - idm: The id of the molecule that was kept (first element of idm_synonyms).

        :param df: The input DataFrame.
        :type df: pandas DataFrame
        :return df_synonyms: The output DataFrame containg synonyms for identifying duplicate molecules.
        :rtype: pandas.DataFrame

        .. warning:: make sure there is no duplicate idm before using this function, as the final filtering of duplicates is performed with a whitelist based on idm.

        .. todo:: optimize code (table format for reference file, numpy for setting up lists when grouping?)
        """
        # check on col_mol
        if self.col_mol not in df.columns:
            raise ValueError(f"Error! No column {self.col_mol} found for col_mol parameter.")
        # check on on
        if self.on != 'inchikey':
            raise ValueError(f"Error! Unauthorized value for on parameter ({self.on}).")
        elif self.on not in df.columns:
            logging.warning(f"Column {self.on} not found, so computing it.")
            df = df.copy()  # this certainly removes the pandas warnings but is it the right way of doing this?
            df.loc[:, self.on] = df.loc[:, self.col_mol].map(rdinchi.MolToInchiKey)
        # init
        df.index = df[self.on]
        df.drop(self.on, axis=1, inplace=True)
        # define synonyms in current dataframe
        df_synonyms = pd.DataFrame(df.groupby(self.on)[self.col_id].apply(list))
        df_synonyms.rename({self.col_id: self.col_id_synonyms}, axis=1, inplace=True)
        df_synonyms[self.col_id] = df_synonyms[self.col_id_synonyms].map(lambda x: x[0])

        # df_synonyms[self.col_id] = df_synonyms[self.col_id_synonyms].map(lambda x: x[0])
        # use information stored in ref file as well, if provided
        if self.ref_file is not None:
            key = Path(self.ref_file).stem
            # init ref file if does not exist already
            if not Path(self.ref_file).is_file():
                self.init_ref_file()
            # open it with a lock as we'll need to update it at the end
            with utils.SafeHDF5Store(self.ref_file) as store:
                df_ref = store[key]
                df_ref = pd.concat([df_ref, df_synonyms])
                df_ref = pd.DataFrame(df_ref.groupby(self.on)[self.col_id_synonyms].apply(list))
                df_ref[self.col_id_synonyms] = df_ref[self.col_id_synonyms].map(lambda x: list(chain.from_iterable(x)))
                df_ref[self.col_id] = df_ref[self.col_id_synonyms].map(lambda x: x[0])
                df_ref.to_hdf(self.ref_file, key=key)
            return df_ref
        else:
            return df_synonyms

    def init_ref_file(self):
        """Initiate an empty reference hdf for identifying duplicates.

        :return: True if the reference file could be initialized, False otherwise.
        :rtype: bool
        """
        try:
            # delete file if it already exists
            if Path(self.ref_file).is_file():
                Path.unlink(self.ref_file)
            # init
            key = Path(self.ref_file).stem
            # create a new ref file
            df_ref = pd.DataFrame({self.on: [], self.col_id_synonyms: [], self.col_id: []})
            df_ref.index = df_ref[self.on]
            df_ref.drop(self.on, axis=1, inplace=True)
            df_ref.to_hdf(self.ref_file, key=key)
            return True
        except ValueError:  # certainly not the only kind of error but PEP8 is against using just plain except.
            return False

    def mark_dupl(self, df):
        """Mark duplicate entries found in df or in ref. Designed to work from
        within Standardizer.run only as it updates columns 'status' and 'task' to
        respectively 'filter' and 'filter_duplicates'. Might still be useful if exposed.


        Moreover, update the reference file (if specified) by using a lock (from utils module) so that
        no other processes can modify it at the same time. This enables a safe duplicate
        molecules filtering accross chunks without having to gather everything in a single file.

        :param df: The input DataFrame.
        :type df: pandas.DataFrame
        :return df_marked: The output DataFrame.
        :type df_marked: pandas.DataFrame
        """
        # any molecule whose id is not in df_synonyms[col_id] is a duplicate of another
        df_synonyms = self.find_synonyms(df)
        df = df.copy()  # supress warnings, but might significantly slow down the process..?
        df.loc[~df.loc[:, self.col_id].isin(df_synonyms.loc[:, self.col_id]), ['status', 'task']] = ('filter', 'filter_dupl')
        return df


class Standardizer(Filter):
    """A class for standardizing molecular structures. The standardization itself is based
    on a protocol that the user can modify.

    By default this protocol consists in 12 tasks applied to each molecule invidually:

        1) **sanitize**: structural check to make sure molecules can be used safely
        2) **disconnect_metal**: break bonds involving metallic atoms, resulting in potentially several molecules per structure.
        3) **keep_largest**: retrieve only the largest molecule (molecular weight in Da) from a structure. Medchem elements are preferred over others, independant from size.
        4) **filter_hac**: filter molecules with a heavy atom count not in the accepted range. By default: hac > 3.
        5) **filter_molweight**: filter molecules with a molecular weight not in the accepted range. By default: molweight <= 1000.0.
        6) **filter_nrings**: filter molecules with a number of rings (Smallest Sets of Smallest Rings or SSSR) not in the accepted range. By default: nrings > 0.
        7) **filter_medchem**: filter molecules with elements not considered as medchem. By default: elements in H, B, C, N, O, F, P, S, Cl, Br, I.
        8) **remove_isotopes**: set all atoms to their most common isotope (i.e. 14C becomes 12C which is C).
        9) **normalize**: always write the same functional groups in the same manner.
        10) **uncharge**: remove all charges on a molecule when it is possible. This is different from rdkit.Chem.MolStandardize.charge module as there is no attempt for reaching the zwitterion.
        11) **canonicalize**: enumerate the canonical tautomer.
        12) **remove_stereo**: remove all remaining stereochemistry flags on the molecule.

    Finally, duplicate entries can be filtered using InChiKeys when postprocessing the DataFrame with molecules (argument in the run_df method).

    This results in new columns in the input DataFrame:

        - the 'mol' column: updated structure (only for the protocol)
        - the 'status' column: either standardized, filter or error.
        - the 'task' column: the latest task that was applied to the molecule.

    The global variable DESCRIPTORS contains keywords for computing descriptors used for filters.
    Although this could be expanded quite easily if needed, there is currently only four descriptors implemented:

        - molweight: molecular weight (Da)
        - hac: heavy atom count
        - nrings: number of rings (Smallest Sets of Smallest Rings or SSSR)
        - elements: atomic symbols

    Two different kind of operations can be defined:

        - numeric: use of comparison operators such as <, >, <=, >=, == and != for numeric comparisons.
        - set: use of ' in ' or ' not in ' (spaces are required!) for inclusions/exclusions.

    Examples of possible operations:

        >>> '0 < molweight <= 1000.0'
        >>> 'hac > 3'
        >>> 'nrings != 0'
        >>> 'elements in H,C,N,O'


    .. note:: For now, the user can only change the task order and the values of filters, but it would be relatively easy to add more functionality.

    .. todo:: Add an argument to specify a configuration file in json.

    .. todo:: Rename standardized status to ready or passed, something shorter. Also filter should be filtered.

    .. todo:: Check latest publication sent by Prof. H. Waldmann, there might be an open-source tool for deglycosylating structures, which could become a new task.
    """

    def __init__(self,
                 protocol=None,
                 col_mol='mol',
                 col_id='idm',
                 filter_duplicates=True,
                 ref_file=None,
                 on='inchikey',
                 ref_dataset=None,
                 suffix=None,
                 elements_medchem={'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'}
                 ):
        """Create a Standardizer object."""
        # filter
        super(Standardizer, self).__init__()
        # standardizer
        self.elements_medchem = elements_medchem
        if protocol is None:
            self._protocol = {'tasks': ['sanitize',
                                        'disconnect_metal',
                                        'keep_largest',
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
        # workers
        self.metal_disconnector = MetalDisconnector()
        self.normalizer = Normalizer()
        self.full_uncharger = FullUncharger()
        self.canonicalizer = TautomerCanonicalizer()

    @property
    def protocol(self):
        return self._protocol

    @protocol.setter
    def protocol(self, protocol):
        # input is a json file => convert it to a dict
        if isinstance(protocol, str):
            utils.check_arg_input_file(protocol)
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

    def remove_isotopes(self, mol):
        """Return a molecule without any isotopes.

        :param mol: The input molecule.
        :type mol: : rdkit.Chem.Mol
        :return mol: The molecule without isotope.
        :rtype mol: : rdkit.Chem.Mol
        """
        for a in mol.GetAtoms():
            a.SetIsotope(0)
        return mol

    def keep_largest(self, mol):
        """Return the largest molecule found in a molecular structure. Molecules with
        only elements usually used in medicinal chemistry are preferred over others,
        no matter their molecular weight.
        The list of medchem elements is contained in the protocol paramater of the Standardizer object:
        (protocol['filter_medchem']).

        :param mol: The input molecule(s).
        :type mol: rdkit.Chem.Mol
        :param store_molweight: Store 'molweight' property in molecule (as str, with precision=4) so it can be used later.
        :type store_molweight: bool
        :return mol: The largest molecule of the input structure.
        :rtype mol: rdkit.Chem.Mol
        """
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        # no need to look further if we have only one fragment!
        if len(frags) < 2:
            return mol
        # otherwise, we have to compare the fragments
        # init
        largest_molweight = -1.0  # so we are sure to update this on the first iteration
        largest_frag = None
        largest_is_medchem = False
        # begin
        for frag in frags:
            # is_medchem
            is_medchem = self.filter_mol(frag, f'elements in {", ".join(str(x) for x in self.elements_medchem)}')
            # molweight
            molweight = Descriptors.ExactMolWt(frag)
            # compare to the current largest fragment
            update_largest = False
            if not largest_is_medchem:
                if is_medchem:
                    update_largest = True
                else:
                    if molweight > largest_molweight:
                        update_largest = True
            else:
                if not is_medchem:
                    pass
                else:
                    if molweight > largest_molweight:
                        update_largest = True
            # update largest
            if update_largest:
                largest_is_medchem = is_medchem
                largest_frag = frag
                largest_molweight = molweight
        # store molweight in case we need it later
        return largest_frag

    @timeout_decorator.timeout(TIMEOUT)
    def _run(self, mol):
        """Helper function for run.
        Contains all tasks defined within the protocol. Since some operations are
        expansive and could last a very long time for complex molecules (normalize, canonicalize),
        a timeout value is set globally. The run function is the one that can catch the exception
        raised by timeouts.

        :param mol: The input molecule.
        :type mol: rdkit.Chem.Mol
        :return result: A tuple containing the molecule, its status and the further task name it reached.
        :rtype: tuple
        """

        # initiate_mol
        if mol is None:
            return (mol, 'error', 'initiate_mol')

        # begin protocol
        mol = copy.deepcopy(mol)
        for task in self._protocol['tasks']:
            # metals
            if task == 'disconnect_metal':
                try:
                    mol = self.metal_disconnector.disconnect(mol)
                except ValueError:
                    return (mol, 'error', 'disconnect_metal')

            # fragments
            elif task == 'keep_largest':
                try:
                    mol = self.keep_largest(mol)
                except ValueError:
                    return (mol, 'error', 'keep_largest')

            # filters
            elif task == 'filter_medchem' or task == 'filter_hac' or task == 'filter_molweight' or task == 'filter_nrings':
                try:
                    if not self.filter_mol(mol, self._protocol[task]):
                        return (mol, 'filter', task)
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
        return (mol, 'standardized', 'passed')

    def run(self, mol):
        """Execute the standardization protocol on a molecule.
        Molecule that exceed the timeout value are filtered with a task='timeout'.

        :param mol: The input molecule.
        :type mol: rdkit.Chem.Mol
        :return result: A tuple containing the molecule, its status and the further task name it reached.
        :rtype: tuple
        """
        try:
            return self._run(mol)
        except timeout_decorator.TimeoutError:
            return (mol, 'filter', 'timeout')

    def run_df(self, df):
        """Apply the standardization protocol on a DataFrame, with the possibility of directly filtering duplicate entries as well.
        This can be very useful as the standardization process can expose duplicate entries due to salts removal, neutralization,
        canonical tautomer enumeration, and stereochemistry centers unlabelling

        If a reference file is specified, duplicate removals becomes possible accross chunks.


        :param mol: The input molecule.
        :type mol: rdkit.Chem.Mol
        :return results: Three DataFrames separated by status:

            - standardized
            - filtered
            - error

        :rtype: tuple

        .. note:: As a side effect, the output DataFrames get indexed by idm and standardized entries keep their 'inchikey' column (mandatory for now).
        """
        # run standardization protocol
        df.index = df[self.col_id]
        df[self.col_mol], df['status'], df['task'] = zip(*df[self.col_mol].map(self.run))
        # do not apply filter duplicates on molecules with errors or that were already filtered for x reasons
        df_error = df[df['status'] == 'error']
        df_filtered = df[df['status'] == 'filter']
        df = df[df['status'] == 'standardized']
        # filter duplicates
        if self.filter_duplicates:
            dupl_filter = DuplicateFilter(on=self.on, col_mol=self.col_mol, col_id=self.col_id, ref_file=self.ref_file)
            df = dupl_filter.mark_dupl(df)
            # clean up for postprocess
            df[dupl_filter.on] = df.index
            df.index = df[self.col_id]
            df.drop(self.col_id, axis=1, inplace=True)
            # separate dupl from the rest
            df_filtered = pd.concat([df_filtered, df[df['status'] == 'filter']], join='inner')  # drop inchikey col as the info is stored in ref_file anyway
            df = df[df['status'] == 'standardized']
        # tuple of dataframes
        return (df, df_filtered, df_error)
