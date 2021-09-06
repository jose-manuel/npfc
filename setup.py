import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()


# begin

setuptools.setup(
    name="npfc",
    version='0.7.13',
    author="Jose-Manuel Gally",
    author_email="josemanuel.gally@mpi-dortmund.mpg.de",
    description="A tool for describing Natural Product- (NP) fragments combinations and identifying pseudo-NPs.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/jose-manuel/npfc",
    packages=setuptools.find_packages(),
    install_requires=[
    'adjustText >= 0.7.3',
    'networkx >= 2.6',
    'ipython >= 7.25',
    'more-itertools >= 8.8',
    'pandas >= 1.1.4',
    'pillow >= 8.3.1',
    'psycopg2-binary >= 2.9',
    #'tables >= 3.6',  # version is older than in conda, and during testing this throws: ValueError: The file 'tests/tmp/test_dupl_ref.hdf' is already opened.  Please close it before reopening.  HDF5 v.1.8.5-patch1, FILE_OPEN_POLICY = 'strict'
    'rdkit-pypi >= 2021.03',
    'snakemake >= 5.0',
    'seaborn >= 0.11',
    'toolz >= 0.11',
    'scipy >= 1.5',
    # extras
    'isort >= 4.3',
    'pylint >= 2.4',
    'pytest >= 6.2',
    'sphinx >= 3.0',
    'sphinx-autodoc-typehints >= 1.10',
    'sphinx_rtd_theme >= 0.4',
    'twine >= 3.4',
    'pytest >= 6.2',
    'wheel >= 0.36',
    ],
    ## the extra requirements below were installed when typing: pip install -e .[dev]
    # extra_require={"dev": ['isort >= 4.3',
    #                        'pylint >= 2.4',
    #                        'pytest >= 6.2',
    #                        'sphinx >= 3.0',
    #                         'sphinx-autodoc-typehints >= 1.10',
    #                         'sphinx_rtd_theme >= 0.4',
    #                         'twine >= 3.4',
    #                         'pytest >= 6.2',
    #                         'wheel >= 0.36',
    #                         ]
    # },
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    keywords=['chemical biology', 'pseudo-natural products', 'computational chemistry', 'natural products', 'fragments'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=[
             # activities
             'bin/act_annotate_fc',
             'bin/act_preprocess',
             'bin/act_preprocess_cpa',
             # chunks
             'bin/chunk_check',
             'bin/chunk_sdf',
             'bin/concat_sdf',
             'bin/concat_synonyms',
             # fc fragments
             'bin/fc/frags_annotate_fcp',
             'bin/fc/frags_search',
              # fc fragment hits
             'bin/fc/fs_filter_fhits',  # rename into fhits_filter
             # fc fragment combinations
             'bin/fc/fc_classify',
             # fc fragment combination graphs
             'bin/fc/fcg_generate',
             'bin/fc/fcg_annotate_pnp',
             # fc molecules
             'bin/fc/mols_concat',
             'bin/fc/mols_dedupl',
             'bin/fc/mols_draw',
             'bin/fc/mols_depict',
             'bin/fc/mols_extract_murcko',
             'bin/fc/mols_load',
             'bin/fc/mols_standardize',
             'bin/fc/mols_subset',
             # fc report
              'bin/fc/report/mols_count',
              'bin/fc/report/report_time',
             # fc commands
             'bin/report_protocol',
             'bin/fc/run_protocol_fc',
             # refs
             'bin/refs_group',
             # helpers
             'bin/peek_hdf',
             ],
    # package_data={'npfc': ['data/mols_deglyco.knwf']},  # key: package name, value: list of data files
    include_package_data=True,
)
