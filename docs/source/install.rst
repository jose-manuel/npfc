============
Installation
============

Currently, the package can be either installed by cloning the github repository
or via PyPi (recommanded).

PyPi
****

In your Python environment, run:

>>> pip install npfc

.. warning: There is currently an issue when pytables is installed via pip.
    We could not investigate this yet, so it currently is strongly recommended
    to install pytables using conda instead.
    To reflect this, we removed pytables from the PyPi package dependecies.
    If not installed, the deduplication step using a reference file will fail.

As a work-around for now, we suggest to create a conda environment and install
the npfc package inside of it:

>>> conda create -n npfc_env python pytables graphviz
>>> conda activate npfc_env
>>> pip install npfc
