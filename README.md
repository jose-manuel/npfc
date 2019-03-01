# npfc

npfc stands for natural product fragment combinations. It is used for classifying
those into predefined categories. It is composed of a package with classes and
tests as well as a main to run the workflow.

## Workflow

The npfc workflow can be broken down in 6 successive steps:
    1) convert a molecular file into a hdf file(s)
    2) standardize structures within the hdf file(s)
    3) remove duplicate entries
    4) search the molecules for requested fragments
    5) classify fragments combinations
    6) generate a report

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


## Stuff to think about

- add extlinks to RDKit documentation (~ MolVS)
- tweak Sanitizeflags to mimic KNIME Mol2RDKit node behavior (partial sanitization).
- implement a git hook so that only changes that pass all tests are committed.
- implement useful logs with at the module level
- build up a workflow using snakemake
- add more checks for inputs for classes/functions
- automatic __version__ var from git tag (versioneer)
- harmonize git commits (commitizen)
- automatic release changelog (standard-version)
- install and use sphinx-autodoc-typehints for much nice docstrings
