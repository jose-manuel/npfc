# npfc

npfc stands for Natural Product Fragment Combinations. It is used for classifying
those into predefined categories. It is composed of a package with classes and
tests as well as a main to run the workflow.

## Workflow

The npfc workflow can be broken down in 6 successive steps:
    1) load/chunk an input file (SDF)
    2) deglycosylate molecules using a KNIME workflow
    3) standardize/filter structures, remove duplicates
    4) for chembl: filter out natural products from chembl (chembl_synth)
    5) search the molecules for requested fragments
    6) classify fragments combinations
    7) create molecular fragment networks

## How to use this

### Inputs

The user requires two main inputs for the workflow:
    1) a molecular file with substructures to use (in this scope: natural products
       or natural product fragments)
    2) a molecular file with molecules to look substructures in

### Outputs

The npfc workflow produces outputs after every task, either in hdf, csv, log or
pdf formats. Please have a look at each step specifically for more details.


## References used for packaging this project

### Package syntax

https://packaging.python.org/tutorials/packaging-projects/
https://guides.github.com/features/mastering-markdown/
https://www.linode.com/docs/applications/project-management/how-to-create-a-private-python-package-repository/

### Automated tests

https://docs.pytest.org/en/latest/goodpractices.html
https://github.com/pluralsight/intro-to-pytest/blob/master/tests/02_special_assertions_test.py
https://github.com/pre-commit/pre-commit/issues/761
https://longair.net/blog/2011/04/09/missing-git-hooks-documentation/
https://pytest.readthedocs.io/en/2.8.7/skipping.html

### Sphinx documentation

https://www.pythonforthelab.com/blog/documenting-with-sphinx-and-readthedocs/
https://protips.readthedocs.io/git-tag-version.html
https://github.com/agronholm/sphinx-autodoc-typehints
https://docs.python.org/3/library/typing.html
https://stackoverflow.com/questions/45957615/check-a-variable-against-union-type-at-runtime-in-python-3-6

### Useful links

Pandas warnings
https://stackoverflow.com/questions/42101382/pandas-dataframe-assign-arguments

HDF warnings
https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Pahl_NotebookTools_Tutorial.ipynb
http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html

Optimizations
https://stackoverflow.com/questions/22219004/grouping-rows-in-list-in-pandas-groupby

Better logging
https://stackoverflow.com/questions/10973362/python-logging-function-name-file-name-line-number-using-a-single-file

Changelog
https://stackoverflow.com/questions/3523534/good-ways-to-manage-a-changelog-using-git

Postgresql queries optimization
https://www.postgresql.org/docs/8.3/queries-limit.html

## Useful tips for building this repo

From the root of the repo (src), type:
python setup.py test

This command actually runs all the tests and execute the following commands:
    1) Create an archive from the package
    python setup.py sdist

    2) Install the package
    cd dist
    pip install npfc-gally-0.0.1.tar.gz

    3) Uninstall the package
    pip uninstall npfc-gally

## Testing the package while developing

Go to src directory and type:
python -m pytest tests -svv

## Generating sphinx docs

Go to src/docs directory and type:
make html
Then open the file docs/build/index.rst with firefox and enjoy!

## Manage package version

Use the command commit found in bin folder. If commits were pushed without any
bump in version, then a suffix is added to the version number (i.e.e 0.0.20-4-g1ac82ee).
The setup.py extracts the version tag from the git repo and the init.py defines the
package version using the setup.py file.

The script commit provides an useful interface for committing changes as it automates several tasks:

    1) add changes for commit
    2) use commitizen for formatting consistent commit messages
    3) define a version name based on its previous value. If upgrade, then the new version
       will overwrite the old version. If not, a "dirty" suffix will be appended to the
       current version.
    4) uninstall a previous installation of npfc, if any
    5) create an archive tar.gz for the package npfc.
    6) install it locally
    7) send the archive to the cluster.

Also, the atom plugin git-log can be summoned with ctrl+shift+p to visualize the history of versions.

## Snakemake

To create the task tree:
>>> snakemake --dag | dot -Tsvg > dag.svg


## Stuff to think about

- add extlinks to RDKit documentation (~ MolVS)
- tweak Sanitizeflags to mimic KNIME Mol2RDKit node behavior (partial sanitization) or split molblocks/smiles before conversion to RDKit
- implement a git hook so that only changes that pass all tests are committed
- automatic release changelog (standard-version?)
- ref_file for duplicate entries removal should be defined in the protocol
- check if ref_file would be better in table mode for appending data instead of rewriting the whole file for each chunk

## Some weird errors encountered during development
    - cannot set user-defined timeout for standardization because the timeout value is set during the loading of the library
    - deglycosylation: KNIME CDK node Sugar Remover is very good at surprising me:
        - loses stereochemistry information
        - if output is smiles, remove all hydrogens so heavy atoms are radicals
        - high failure rate (10%) which returns empty molecules
        - if the sugar is inside of the molecule, then fragments are returned. Yet to define if the larger fragment is the better one
