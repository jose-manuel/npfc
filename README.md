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

### Conda environment

Export conda env to yml file:

>>> conda env export > conda_npfc.yml

Use this yml file to create a conda env:

>>> conda env create -f conda_npfc.yml

Careful as two dependencies were installed with pip and require special handling:

#### timeout_decorator

This package is used for setting a timeout for standardizing molecules.

To install it, just run:

>>> pip install timeout_decorator

#### pdbeccdutils

This package is used for comparing 2D depictions of a same molecule and retrieve
the "best" one.

Unfortunately, the latest release (0.5, tested on the 28th of May 2019) is broken and
cannot be imported anymore in my environment (missing "dataclasses"). So I just used
the 0.4 git clone I have on my laptop instead. I sent this copy to the cluster.
This is very whacky and I hope I can find a better way before releasing this package.

Inside of the package root folder:

>>> python setup.py sdist
>>> cd dist
>>> pip install pdbeccdutils-0.4.tar.gz

Or much more easier:

>>> pip install git+https://gitlab.ebi.ac.uk/pdbe/ccdutils.git

#### networkx

I don't know why networkx is a special case. It is installed using conda install.

Since I intend to have the same version everywhere:

>>> conda install networkx==2.2


#### knime

Things got out of hand and I now have to run a KNIME worklow in the pipeline.


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
https://stackoverflow.com/questions/25866102/how-do-we-embed-images-in-sphinx-docs
https://stackoverflow.com/questions/22149669/how-can-i-place-images-side-by-side-in-restructured-text

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

>>> python -m pytest tests -svv

Also, to skip local workflow tests:

>>> python -m pytest --ignore-glob='*06*' tests -sv

## Generating sphinx docs

Go to src/docs directory and type:
make html
Then open the file docs/build/index.rst with firefox and enjoy!

## Generating the class hierarchy diagram

From the project root folder, type:

(npfc) gally@m18047-lin:~/Projects/NPFC/src$ pyreverse -my -A -o svg -p npfc npfc/

This creates a file named classes_npfc.svg.


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

### Useful Commands

To create the task tree:
>>> snakemake --dag | dot -Tsvg > dag.svg

To run a Snakemake job locally:
>>> snakemake -s run_crms.smk

To run a Snakemake on the cluster:
>>> snakemake -s run_crms.smk --jobs 999 --cluster 'bsub -q mpi-short -R "scratch" -J {cluster.name}'

### NPFC pipelines

Actually, 3 pipelines are present in the NPFC project:
    - prep_crms:
    - fcc_dnp
    - fcc_chembl


## Run jobs on clusters

### GWDG1

The Goettingen cluster (Gesellschaft fuer wissenschaftliche Datenverarbeitung mit mbH Goettingen)
can be used for running jobs using either lsf (bsub) or slurm (sbatch).
There is something wrong with the way lsf is configured as job parameters specified in bsub scripts
are simply ignored. This makes that one has to specify everything with command line arguments.


### CLEM

clem is the name of the local slurm cluster at the MPI-Dortmund. It is mainly used
for processing microscope imaging but the nodes e045 to e047 (3 nodes) are Axel's to use for
Cell Painting. When there is no CP going on, I can use the nodes for NPFC.
In order to avoid conflicts with other research groups, I have to specifiy the correct
queue (--partition comas). Also if I want to run one job per core (24 cores per node),
I also have to use the argument --oversubscribe.

It is also possible to check how many cores are being used by directly logging into the
node from the front-end and use htop.

A ganglia server is also accessible from the browser:

http://clem.mpi-dortmund.mpg.de/ganglia/


## Stuff to think about

- add extlinks to RDKit documentation (~ MolVS)
- tweak Sanitizeflags to mimic KNIME Mol2RDKit node behavior (partial sanitization) or split molblocks/smiles before conversion to RDKit
- implement a git hook so that only changes that pass all tests are committed
- automatic release changelog (standard-version?)
- ref_file for duplicate entries removal should be defined in the protocol
- check if ref_file would be better in table mode for appending data instead of rewriting the whole file for each chunk
- problem with the way local workflows are tested: cannot be tested on the fly as they use the currently installed npfc version

## Some weird errors encountered during development
    - cannot set user-defined timeout for standardization because the timeout value is set during the loading of the library
    - deglycosylation: KNIME CDK node Sugar Remover is very good at surprising me:
        - loses stereochemistry information
        - if output is smiles, remove all hydrogens so heavy atoms are radicals
        - high failure rate (10%) which returns empty molecules
        - if the sugar is inside of the molecule, then fragments are returned. Yet to define if the larger fragment is the better one

## TODO

    - fragment: count the number of combinations involving repeating fragments, store it into graph as weight for edges, only connection type will be used for pnp annotation
    - automate the nfn computation and display using cytoscape.js/CyREST?
    - automate the rendering of MFGs using cytoscape.js/CyREST?
    - annotate molecules with bioactivities
    - standardize: switch to csv file and add rotate ref file using an option for max num of entries (1M?)
