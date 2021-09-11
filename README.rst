npfc: Natural Product Fragment Combinations
===========================================

**npfc** is a chemoinformatics tool for classifying Natural Product (NP) fragment
combinations into predefined categories and therefore identifying pseudo-NPs.

Pseudo-NPs are novel NP-inspired compound classes that combine the biological
relevance of NPs with the efficient exploration of chemical space by
fragment-based drug design.

The npfc tool is written in Python and based on several key packages:

- `RDKit`_ for handling chemistry
- `pandas`_ for managing data into DataFrames
- `NetworkX`_ for modelling graphs
- `Snakemake`_ for encapsuling scripts into reproducible workflows

Installation
============

The npfc tool can be installed using PyPi. In your Python environment, run:

>>> pip install npfc

Warning!!! There is currently an issue when pytables is installed via pip.
It was removed from the dependency but is required for the npfc workflow
(reference file used during deduplication). Instead, please install pytables via
conda:

>>> conda create -n npfc_env python pytables
>>> conda activate npfc_env
>>> pip install npfc

Documentation
=============

The full documentation is available at: https://github.com/jose-manuel/npfc.
It describes the API, as well as the different workflows implemented.

Contribution
============

Feedback from the community is warmly welcomed. It can be in the form of bug
reports and feature requests submitted via github or code contribution via
forking this repo and submitting pull requests.

License
=======

npfc is licensed under the `MIT license`_.

.. _`RDKit`: http://www.rdkit.org
.. _`pandas`: https://pandas.pydata.org/
.. _`NetworkX`: https://networkx.org/
.. _`Snakemake`: https://snakemake.readthedocs.io/en/stable/
.. _`MIT license`: https://github.com/jose-manuel/npfc/blob/master/LICENSE
