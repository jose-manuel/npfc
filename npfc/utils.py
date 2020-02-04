"""
Module utils
================
"""

# standard
import logging
from pathlib import Path
# data handling
import pickle
import base64
# docstrings
from rdkit.Chem import Mol
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdinchi import MolToInchiKey
from typing import Union
from typing import List

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# allowed suffixes
EXTS_INPUT = [['.sdf'], ['.sdf', '.gz'],
              ['.csv'], ['.csv', '.gz'],
              ['.hdf']]  # , ['.feather']]  # tests with feather work but not in production, dig into that later

EXTS_CONFIG = [['.json']]

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
    logging.debug(f"Input file Suffixes are: {suffixes}")
    # is the file an archive?
    if len(suffixes) > 1:
        compression = suffixes[1]
        logging.debug(f"File is an archive with compression={compression}")
        # special case for gzip so .gz files can be read directly with pandas
        if compression == '.gz':
            compression = 'gzip'
        else:
            raise ValueError(f"Error! Unexpected value for compression suffix: '{compression}'.")
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


def check_arg_positive_number(value: Number) -> bool:
    """Return True of the value is indeed a positive number (>0), raise a TypeError otherwise.

    :param value: the argument to test
    """
    logging.debug(f"value is {value} ({type(value)})")
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


def _configure_logger(log_level: str) -> logging:
    """Configure the logging in a centralized way. This is useful for scripts mostly.

    :param log_level: the logging level to use, accepted values are: CRITICAL, ERROR, WARNING, INFO, DEBUG.
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')
    logging.basicConfig(level=numeric_level, format='%(asctime)s -- %(levelname)s -- %(message)s')
    logger = logging.getLogger(__name__)
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
