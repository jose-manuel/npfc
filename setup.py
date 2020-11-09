import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


# begin

setuptools.setup(
    name="npfc",
    version='0.7.0-6-geb4b3a0',
    author="Jose-Manuel Gally",
    author_email="josemanuel.gally@mpi-dortmund.mpg.de",
    description="A package for analyzing fragment combinations in natural and synthetic molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=[
          'pip',
      ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    keywords=['chemical biology', 'pseudo-natura products', 'chemoinformatics', 'natural products', 'fragments'],
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
             # fct nodes
             'bin/fct_assay',
             'bin/fct_dataset',
             'bin/fct_fragment',
             'bin/fct_molecule',
             # fct relationships
             'bin/fct_fragment_fragment',
             'bin/fct_fragment_molecule',
             'bin/fct_molecule_dataset',
             'bin/fct_molecule_molecule',
             # fct commands
             'bin/fct_run_protocol',
             # fragments
             'bin/frags_annotate_fcp',
             'bin/frags_search',
             # fragment hits
             'bin/fs_filter_fhits',  # rename into fhits_filter
             # fragment combinations
             'bin/fc_classify',
             # fragment combination graphs
             'bin/fcg_generate',
             'bin/fcg_annotate_pnp',
             # molecules
             'bin/mols_concat',
             'bin/mols_dedupl',
             'bin/mols_draw',
             'bin/mols_depict',
             'bin/mols_extract_murcko',
             'bin/mols_load',
             'bin/mols_standardize',
             'bin/mols_subset',
             # protocols
             'bin/report_protocol',
             'bin/run_protocol',
             # helpers
             'bin/peek_hdf',
             ],
    # package_data={'npfc': ['data/mols_deglyco.knwf']},  # key: package name, value: list of data files
    include_package_data=True,
)
