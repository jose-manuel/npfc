"""
Module utils
================
"""

# standard
import logging
from pathlib import Path
import os
import time
# data science
from pandas import HDFStore
# docstrings
from typing import Union
from typing import List

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# allowed suffixes
EXTS_INPUT = [['.sdf'], ['.sdf', '.gz'], ['.sdf', '.zip'],
              ['.csv'], ['.csv', '.gz'], ['.csv', '.zip'],
              ['.hdf']]

EXTS_CONFIG = [['.json']]

# types
Number = Union[int, float]
Output_files = List[List[Union[str, int]]]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


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


def check_arg_input_file(input_file: str) -> bool:
    """Return True of the input_file exists, raise an error otherwise.

    :param input_file: the input file
    :param input_format: the expected format of the input file
    """
    path_input_file = Path(input_file)
    if not path_input_file.is_file():
        raise ValueError(f"Error! Input file could not be found at {input_file}.")
    if path_input_file.suffixes not in EXTS_INPUT:
        raise ValueError(f"Error! Unexpected '{path_input_file.suffixes}' for input format.")

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
            logging.warning(f"Output_dir could not be found at {output_dir}, attempting to create it.")
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise ValueError(f"Error! Output_dir could not be found at {output_dir}.")

    return True


def check_arg_config_file(config_file: str) -> bool:
    """Return True of the config_file exists, raise an error otherwise.


    :param input_file: the input file
    :param input_format: the expected format of the input file
    """
    path_config_file = Path(config_file)
    if not path_config_file.is_file():
        raise ValueError(f"Error! Input file could not be found at {config_file}.")
    if path_config_file.suffixes not in EXTS_CONFIG:
        raise ValueError(f"Error! Unexpected '{path_config_file.suffixes}' for config format.")

    return True

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class SafeHDF5Store(HDFStore):
    """Implement safe HDFStore by obtaining file lock. Multiple writes will queue if lock is not obtained.

    Copied from https://stackoverflow.com/questions/41231678/obtaining-a-exclusive-lock-when-writing-to-an-hdf5-file.
    """

    def __init__(self, *args, **kwargs):
        """Initialize and obtain file lock."""
        interval = kwargs.pop('probe_interval', 1)
        self._lock = "%s.lock" % args[0]
        while True:
            try:
                self._flock = os.open(self._lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                break
            except (IOError, OSError):
                time.sleep(interval)

        HDFStore.__init__(self, *args, **kwargs)

    def __exit__(self, *args, **kwargs):
        """Exit and remove file lock."""
        HDFStore.__exit__(self, *args, **kwargs)
        os.close(self._flock)
        os.remove(self._lock)
