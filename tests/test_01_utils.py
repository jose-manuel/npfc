"""
Module test_01_utils
====================
Tests for the utils module.

Uses the files generated by the test_00_save.
"""
# standard
from pathlib import Path
# tests
import pytest
# dev library
from npfc import utils

# import logging
# logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_shs_init():
    """Test if a SafeHDF5Store object can be instanciated and the lock is used."""
    store = 'tests/tmp/store.hdf'
    # remove an eventual lock from previous run
    path_lock = Path(store + '.lock')
    if path_lock.is_file():
        path_lock.unlink()
    # open a store object using with allows for automatic exit. If exit is not performed,
    # then core dumps due to stack overflows can happen. This can be prevented by deleting
    # the lock file manually.
    with utils.SafeHDF5Store(store) as STORE:
        # instance of SafeHDF5Store
        assert isinstance(STORE, utils.SafeHDF5Store) is True
        # lock created
        assert path_lock.is_file() is True
    # lock removed
    assert path_lock.is_file() is False


def test_check_arg_bool():
    """Test the parsing of arguments that should be bool only."""
    # should pass
    assert utils.check_arg_bool(True)
    # should throw a TypeError
    error = False
    try:
        assert utils.check_arg_bool(0)
    except TypeError:
        error = True
    assert error is True


def test_check_arg_positive_number():
    """Test the parsing of arguments that should be positive integers only."""
    # should pass
    assert utils.check_arg_positive_number(1)
    # should throw a ValueError
    error = False
    try:
        assert utils.check_arg_positive_number(0)
    except ValueError:
        error = True
    assert error is True
    # should throw a TypeError
    error = False
    try:
        assert utils.check_arg_positive_number('a')
    except TypeError:
        error = True
    assert error is True


def test_check_arg_input_file():
    """Test the parsing of arguments that should be a input file."""
    # should pass
    utils.check_arg_input_file('tests/tmp/test_save_simple.csv')
    utils.check_arg_input_file('tests/tmp/test_save_simple.hdf')
    utils.check_arg_input_file('tests/tmp/test_save_simple.sdf')
    utils.check_arg_input_file('tests/tmp/test_save.sdf.gz')
    utils.check_arg_input_file('tests/tmp/test_save.csv.zip')
    # should throw a ValueError
    error = False
    try:
        utils.check_arg_input_file('tests/tmp/test_file_does_not_exist.csv')
    except ValueError:
        error = True
    assert error is True


def test_arg_output_file():
    """Test the parsing of arguments that should be an output file."""
    # should pass
    utils.check_arg_output_file('tests/tmp/test_new_file.csv')
    utils.check_arg_output_file('tests/tmp/test_new_file.csv.zip')
    utils.check_arg_output_file('tests/tmp/test_new_file.hdf')
    utils.check_arg_output_file('tests/tmp/test_new_file.sdf')
    utils.check_arg_output_file('tests/tmp/test_new_file.sdf.gz')
    # should throw a ValueError
    error = False
    try:
        utils.check_arg_output_file('tests/tmp/test_file.csv.other')
    except ValueError:
        error = True
    assert error is True


def test_arg_config_file():
    """Test the parsing of arguments that should be a config file."""
    # should pass
    utils.check_arg_config_file('tests/tmp/std_protocol.json')
    # should throw a ValueError
    error = False
    try:
        utils.check_arg_config_file('tests/tmp/test.csv')
    except ValueError:
        error = True
    assert error is True
