"""
Module utils
================
"""

from typing import Union

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
    if not isinstance(value, Number):
        raise TypeError(f"Error! Expected a positive number but got {type(value)} instead ({value}).")
    elif value <= 0:
        raise ValueError(f"Error! Expected a positive number but got {value}) instead.")
    return True
