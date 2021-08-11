"""
Module utils
================
"""

# standard
import logging
from pathlib import Path
import sys
import itertools
import signal
from contextlib import contextmanager
# data handling
import pickle
import base64
# typing
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Mol
# chemoinformatics
from typing import Union
from typing import List


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# allowed suffixes
EXTS_INPUT = [['.sdf'], ['.sdf', '.gz'],
              ['.csv'], ['.csv', '.gz'],
              ['.log'],  # for loading formatted log files
              ['.hdf'],
              ['.pkl', '.gz'], ['.model', '.gz']  # model files for computing sa and np scores
              ]  # , ['.feather']]  # tests with feather work but not in production, dig into that later

EXTS_CONFIG = [['.json'], ['.mplstyle']]

# types
Number = Union[int, float]
Output_files = List[List[Union[str, int]]]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def get_file_format(input_file: str) -> tuple:
    """Deduce how the file should be parsed based on its suffixes.

    Only compressions for one single file with initial extension still visibile
    work (see example below).

    >>> from pathlib import Path
    >>> from npfc import utils
    >>> utils.get_file_format(Path('file.csv.gz').suffixes)
    >>> # returns ('CSV', 'gzip')
    >>> utils.get_file_format(Path('file.sdf').suffixes)
    >>> # returns ('SDF', None)

    :param suffixes: suffixes of a file
    :return: a tuple with syntax (format, compression)
    """
    suffixes = Path(input_file).suffixes
    logging.debug("Input file Suffixes are: %s", suffixes)
    # is the file an archive?
    if len(suffixes) > 1:
        compression = suffixes[1]
        logging.debug("File is an archive with compression=%s", compression)
        # special case for gzip so .gz files can be read directly with pandas
        if compression == '.gz':
            compression = 'gzip'
        else:
            raise ValueError("Error! Unexpected value for compression suffix: '%s'.", compression)
    else:
        compression = None

    # return format and compression
    if suffixes[0] in ('.sdf', '.mol', '.sd'):
        return ('SDF', compression)
    elif suffixes[0] == '.csv':
        return ('CSV', compression)
    elif suffixes[0] in ('.hdf', 'hf5'):
        return ('HDF', compression)
    elif suffixes[0] == ".feather":
        return ('FEATHER', compression)


def check_arg_bool(value: bool) -> bool:
    """Return True of the value is indeed a boolean, raise a TypeError otherwise.

    :param value: the argument to test
    """
    if not isinstance(value, bool):
        raise TypeError(f"Error! Expected a boolean but got {type(value)} instead ({value}).")

    return True


def parse_argparse_boolstring(value: str) -> bool:
    """Return True or False given the provided string. If the string is actually not a boolean, raise a TypeError.

    :param value: the argument to test
    :return: the parsed boolean as type bool
    """
    if value.lower() in ("yes", "true", "t", "1"):
        return True
    elif value.lower() in ("no", "false", "f", "0"):
        return False
    else:
        raise TypeError(f"Error! Could not parse provided boolean string ({value}).")


def check_arg_positive_number(value: Number) -> bool:
    """Return True of the value is indeed a positive number (>0), raise a TypeError otherwise.

    :param value: the argument to test
    """
    logging.debug("value is %s (%s)", value, type(value))
    # possible value for Number is None if argument is left unset
    if value is None:
        return True
    # if not None, then looks what it is
    if not isinstance(value, Number.__args__):  # might be a hack but Union object is not compatible with isinstance
        raise TypeError(f"Error! Expected a positive number but got {type(value)} instead ({value}).")
    elif value <= 0:
        raise ValueError(f"Error! Expected a positive number but got {value} instead.")

    return True


def check_arg_input_file(input_file: str) -> bool:
    """Return True of the input_file exists, raise an error otherwise.

    :param input_file: the input file
    :param input_format: the expected format of the input file
    """
    path_input_file = Path(input_file)
    if not path_input_file.is_file():
        raise ValueError(f"Error! Input file could not be found at '{input_file}'.")
    if path_input_file.suffixes not in EXTS_INPUT:
        raise ValueError(f"Error! Unexpected '{path_input_file.suffixes}' for input format.")

    return True


def check_arg_input_dir(input_dir: str) -> bool:
    """Return True of the input_dir exists.

    :param input_dir: the output directory
    :return: True if the directory exists.
    """
    # output_format
    p = Path(input_dir)
    if not p.is_dir():
        raise ValueError(f"Error! Input dir could not be found at '{input_dir}'.")

    return True


def check_arg_output_file(output_file: str, create_parent_dir: bool = True) -> bool:
    """Return True of the output_file has the expected format (deduced from the file extension).

    If the parent directory of the output file does not exist, it has to either be created or fail the check.

    :param output_file: the output file
    :param create_parent_dir: create the output file's parent folder in case it does not exist
    """
    # output_format
    path_output_file = Path(output_file)
    if path_output_file.suffixes not in EXTS_INPUT:
        raise ValueError(f"Error! Unexpected value '{path_output_file.suffixes}' for output format.")

    # create_parent_dir
    output_dir = path_output_file.resolve().parent
    if not output_dir.is_dir():
        if create_parent_dir:
            logging.warning(f"Output_dir could not be found at '{output_dir}', attempting to create it.")
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise ValueError(f"Error! Output_dir could not be found at '{output_dir}'.")

    return True


def check_arg_output_plot(output_file: str, create_parent_dir: bool = True) -> bool:
    """Return True of the output_plot has the expected format (deduced from the file extension).

    Accepted extensions are: svg and png.

    If the parent directory of the output file does not exist, it has to either be created or fail the check.

    :param output_file: the output file
    :param create_parent_dir: create the output file's parent folder in case it does not exist
    """
    # output_format
    path_output_file = Path(output_file)
    if path_output_file.suffixes not in [['.svg'], ['.png']]:
        raise ValueError(f"Error! Unexpected value '{path_output_file.suffixes}' for output plot.")

    # create_parent_dir
    output_dir = path_output_file.resolve().parent
    if not output_dir.is_dir():
        if create_parent_dir:
            logging.warning(f"Output_dir could not be found at '{output_dir}', attempting to create it.")
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise ValueError(f"Error! Output_dir could not be found at '{output_dir}'.")

    return True


def check_arg_output_dir(output_dir: str) -> bool:
    """Return True of the output_dir can exist.

    If the parent directory of the output dir does not exist, it has to either be created or fail the check.

    :param output_dir: the output directory
    :param create_parent_dir: create the output file's parent folder in case it does not exist
    """
    # output_format
    path_output_file = Path(output_dir)
    if not path_output_file.is_dir():
        path_output_file.mkdir(parents=True, exist_ok=True)

    return True


def check_arg_config_file(config_file: str) -> bool:
    """Return True of the config_file exists, raise an error otherwise.

    :param input_file: the input file
    """
    if config_file is not None:
        path_config_file = Path(config_file)
        if not path_config_file.is_file():
            raise ValueError(f"Error! Input file could not be found at '{config_file}'.")
        elif path_config_file.suffixes not in EXTS_CONFIG:
            raise ValueError(f"Error! Unexpected '{path_config_file.suffixes}' for config format.")

    return True


def _configure_logger(log_level: str, logger_name: str = None, log_file: str = None, reset_handlers: bool = True) -> logging:
    """Configure the logging in a centralized way. This is useful for scripts mostly.

    :param log_level: the logging level to use, accepted values are: CRITICAL, ERROR, WARNING, INFO, DEBUG.
    :param logger_name: specify a name for the logger. If none, __name__ is used. If within a loop, better speficy a logger name, as logging messages would otherwise stack.
    :param log_file: if specified, the logging messages will be written to the corresponding file, in addition to stdout.
    :return: a logger object that can be reset
    """
    if reset_handlers:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
    if logger_name is None:
        logger_name = __name__
    # define level with DEBUG instead of logging.DEBUG or numbers
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')

    # create the logger
    logger = logging.getLogger(logger_name)
    if len(logger.handlers) > 0:
        logger.handlers.clear()
    logger.setLevel(numeric_level)
    # define format
    if log_level == 'DEBUG':
        format_string = ("%(asctime)s -- %(levelname)s -- %(funcName)s:%(lineno)d -- %(message)s")
    else:
        format_string = ("%(asctime)s -- %(levelname)s -- %(message)s")
    log_format = logging.Formatter(format_string)

    # Creating and adding the console handler
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)

    # Creating and adding the file handler
    file_handler = logging.FileHandler(logger_name, mode='w')
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    # avoid tons of duplicated log messages when executing code from modules
    logger.propagate = False

    return logger


def encode_object(element: object) -> str:
    """Convert an object to a base64 string.

    :param element: an object to encode
    :return: an base64 string upon success, None otherwise
    """
    try:
        return base64.b64encode(pickle.dumps(element)).decode("utf-8")
    except TypeError:
        return None


def decode_object(string: str) -> object:
    """Convert a base64 string to an object.

    :param string: a base64 string encoding an object
    :return: an object upon success, None otherwise
    """
    try:
        return pickle.loads(base64.b64decode(string))
    except TypeError:
        return None


def encode_mol(mol: Mol) -> str:
    """Convert a molecule to a base64 string.

    :param mol: the input molecule
    :return: the molecule in base64
    """
    try:
        return base64.b64encode(mol.ToBinary()).decode("utf-8")
    except AttributeError:
        return None


def decode_mol(string: str) -> Mol:
    """Convert a string to a RDKit Mol object.

    :param string: a string with a Mol object in bytes with a base64 string representation
    :return: a Mol object upon success, None otherwise

    """
    try:
        return Mol(base64.b64decode(string))
    except TypeError:
        return None


def encode_mol_smiles(mol: Mol) -> str:
    """Convert a mol to a Smiles.

    :param string: a string with a Mol object in bytes with a base64 string representation
    :return: a str object upon success, None otherwise

    """
    try:
        return Chem.MolToSmiles(mol)
    except TypeError:
        return None


def decode_mol_smiles(string: str) -> Mol:
    """Convert a Smiles to a Mol.

    :param string: a string with a Mol object in bytes with a base64 string representation
    :return: a Mol object upon success, None otherwise

    """
    if string is None:
        string = ''
    try:
        return Chem.MolFromSmiles(string)
    except TypeError:
        return None


def get_shortest_path_between_frags(mol: Mol, aidxf1: set, aidxf2: set) -> tuple:
    """Return the shortest path within a molecule between two fragments defined by atom indices.
    First and last atom indices are part of respectively fragment 1 and fragment 2, so they should not
    be considered when estimating the distance between fragments.

    (i.e. distance = len(shortest_path) - 2)

    :param mol: The input molecule.
    :param aidxf1: the atom indices of the first fragment found in the molecule
    :param aidxf2: the atom indices of the second fragment found in the molecule
    :return: the atom indices of the shortest path between both fragments. The first index is the attachment point from fragment 1 whereas the last index is the attachment point from fragment 2
    """
    # 1/ compute every pairwise atom combination between both fragments
    pairwise_combinations = itertools.product(tuple(aidxf1), tuple(aidxf2))
    # 2/ for each of those, compute the shortest path possible
    pairwise_combinations = list(pairwise_combinations)
    all_paths = [Chem.GetShortestPath(mol, pc[0], pc[1]) for pc in pairwise_combinations]
    logging.debug("Looking for the shortest path shortest path among these:")
    [logging.debug("Path (%s): %s", str(i).zfill(3), p) for i, p in enumerate(all_paths)]
    # 3/ return one of the shortest pathes
    return min(all_paths, key=lambda x: len(x))


def fuse_rings(rings: tuple) -> list:
    """
    Check for atom indices in common between rings to aggregate them into fused rings.

    :param rings: the ring atoms as provided by the RDKit function mol.GetRingInfo().AtomRings() (iteratble of iteratable of atom indices).
    :return: the fused ring atoms (list of lists of atom indices)
    """
    # condition to exit from recursive fusion of ring atom indices
    done = False

    while not done:
        fused_rings = []
        num_rings = len(rings)
        # pairwise check for common atoms between rings
        for i in range(num_rings):
            # define a core
            fused_ring = set(rings[i])
            for j in range(i+1, num_rings):
                # detect if ring is fused
                if set(rings[i]) & set(rings[j]):
                    # add fused ring to our rign atom list
                    fused_ring = fused_ring.union(rings[j])
            # either lone or fused ring, first check if not already in fused_rings
            if any([fused_ring.issubset(x) for x in fused_rings]):
                continue
            fused_rings.append(list(fused_ring))
        rings = list(fused_rings)
        # there are no rings to fuse anymore
        if num_rings == len(rings):
            done = True

    return rings


@contextmanager
def timeout(time):
    """This function is used within a with statement:

    >>> with timeout(5):
    >>>    do something

    If the code block execution time exceeds the time threshold, a TimeoutError is raised.

    :param time: time in seconds allowd to the code block before aborting its execution

    References

        - https://www.jujens.eu/posts/en/2018/Jun/02/python-timeout-function/
        - https://docs.python.org/3/library/contextlib.html

    """
    # register a function to raise a TimeoutError on the signal.
    signal.signal(signal.SIGALRM, raise_timeout)
    # schedule the signal to be sent after time
    signal.alarm(time)
    # run the code block within the with statement
    try:
        yield
    except TimeoutError:
        pass  # exit the with statement
    finally:
        # unregister the signal so it won't be triggered if the timeout is not reached
        signal.signal(signal.SIGALRM, signal.SIG_IGN)


def raise_timeout(signum, frame):
    """Function to actually raise the TimeoutError when the time has come.
    """
    raise TimeoutError
