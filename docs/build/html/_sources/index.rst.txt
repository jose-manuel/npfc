.. npfc documentation master file, created by
   sphinx-quickstart on Mon Feb 11 10:18:59 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

********************************
Welcome to npfc's documentation!
********************************

This Python 3 package is designed for analyzing Natural-Product Fragment
Combination (npfc).

It is mostly based on the RDKit, Pandas and NetworkX.

On top of the existing modules, a set of scripts are included as well.
The execution of scripts is orchestrated via Snakemake workflows.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   deduplicate
   draw
   filter
   fragment_search
   fragment_combination
   fragment_combination_graph
   fragment_combination_point
   load
   notebook
   save
   standardize
   utils

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

====================
Package Architecture
====================

Module interactions
********************

.. image:: _images/packages_npfc.svg
    :align: center


Class hierarchy
***************

Right-click open in new tab to have a better view of this figure.

.. image:: _images/classes_npfc.svg
    :align: center
