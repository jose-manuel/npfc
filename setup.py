import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


# begin

setuptools.setup(
    name="npfc",
    version='0.5.16-1-g6e86a2a',
    author="Jose-Manuel Gally",
    author_email="josemanuel.gally@mpi-dortmund.mpg.de",
    description="A package for analywing fragment combinations in natural and synthetic molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=[
          'pip',
      ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    keywords=['chemoinformatics', 'natural products', 'fragments', 'chemical biology', 'pseudo-natura products'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=[
             'bin/act_annotate_fc',
             'bin/act_preprocess',
             'bin/act_preprocess_cpa',
             'bin/chunk_sdf',
             'bin/mols_concat',
             'bin/mols_dedupl',
             'bin/mols_deglyco',
             'bin/mols_draw',
             'bin/mols_depict',
             'bin/mols_extract_murcko',
             'bin/mols_load',
             'bin/mols_standardize',
             'bin/mols_subset',
             'bin/mols_substruct',
             'bin/fc_classify',
             'bin/fgraph_generate',
             'bin/fgraph_annotate_pnp',
             'bin/report_protocol',
             'bin/run_protocol',
             'bin/peek_hdf',
             ],
    # package_data={'npfc': ['data/mols_deglyco.knwf']},  # key: package name, value: list of data files
    include_package_data=True,
)
