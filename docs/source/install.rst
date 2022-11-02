============
Installation
============

Currently, the package can be either installed by cloning the github repository
or via PyPi (recommanded).


Conda environment setup
***********************

First, you should setup your conda environment.
One possibility is to create a new conda environment using this :download:`file <_data/npfc_test.yml>`.

>>> conda env create -f npfc_test.yml

Another alternative is the install the dependencies in an existing environment:

>>> conda install python pytables graphviz datrie h5py recommonmark ipykernel


Install the NPFC package
************************

NPFC is available on PyPi.
In your Python environment, run:

>>> pip install npfc

.. warning: There is currently an issue if pytables is installed via pip.
    We could not investigate this yet, so it currently is strongly recommended
    to install pytables using conda instead.
    To reflect this, we removed pytables from the PyPi package dependecies.
    If not installed, the deduplication step using a reference file will fail.

