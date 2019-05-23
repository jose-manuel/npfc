"""
tests.test_06_local_workflow
~~~~~~~~~~
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.
"""
# standard library
import logging
# data handling
import pandas as pd
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc import fragment

# debug
# logging.basicConfig(level=logging.DEBUG)

## TO DO LATER
