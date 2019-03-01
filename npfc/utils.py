"""
Module utils
================
"""

# standard
from pathlib import Path
from typing import Union
from typing import List

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Number = Union[int, float]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
    # possible value for Number is None if argument is left unset
    if value is None:
        return True
    # if not None, then looks what it is
    if not isinstance(value, Number.__args__):  # might be a hack but Union object is not compatible with isinstance
        raise TypeError(f"Error! Expected a positive number but got {type(value)} instead ({value}).")
    elif value <= 0:
        raise ValueError(f"Error! Expected a positive number but got {value} instead.")
    return True


def check_arg_input_file(input_file: str, input_format: List[str] = None) -> bool:
    """Return True of the input_file exists, raise a ValueError otherwise.
    If specified, also checks for the expected format of the file, defined by its extension.
    The extension follows the syntax of the Path().suffixes methods, so 'file.csv' would be ['.csv']
    and 'file.csv.gz' would be ['.csv', '.gz'].

    :param input_file: the input file
    :param input_format: the expected format of the input file
    """
    path_input_file = Path(input_file)
    if not path_input_file.is_file():
        raise ValueError(f"Error! Input file could not be found at {input_file}.")
    if input_format is not None and path_input_file.suffixes != input_format:
        raise ValueError(f"Error! Expected '{input_format}' for input format but got instead {path_input_file.suffixes}.")
    return True
