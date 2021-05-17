import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


# begin

setuptools.setup(
    name="npfc",
    version='0.7.11-2-ge9d8e7d',
    author="Jose-Manuel Gally",
    author_email="josemanuel.gally@mpi-dortmund.mpg.de",
    description="A package for analyzing fragment combinations in natural and synthetic molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=[
          'pip',
        #  'rdkit',
        #  'networkx',
        #  'pandas',
        #  'seaborn',
        #  'toolz',
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
             # fct nodes
             'bin/fct/fct_assay',
             'bin/fct/fct_dataset',
             'bin/fct/fct_document',
             'bin/fct/fct_fragment',
             'bin/fct/fct_molecule',
             'bin/fct/fct_species',
             'bin/fct/fct_target',
             # fct relationships
             'bin/fct/fct_assay_document',
             'bin/fct/fct_assay_species',
             'bin/fct/fct_fragment_fragment',
             'bin/fct/fct_fragment_molecule',
             'bin/fct/fct_molecule_assay',
             'bin/fct/fct_molecule_dataset',
             'bin/fct/fct_molecule_document',
             'bin/fct/fct_molecule_molecule',
             'bin/fct/fct_target_assay',
             'bin/fct/fct_target_species',
             # fct commands
             'bin/fct/report/report_fct_molecule',
             'bin/fct/run_protocol_fct',
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
