import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


# begin

setuptools.setup(
    name="npfc",
    version='0.2.6',
    author="Jose-Manuel Gally",
    author_email="josemanuel.gally@mpi-dortmund.mpg.de",
    description="A package for classifying fragment combinations in molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    keywords=['chemoinformatics', 'natural products', 'fragments', 'chemical biology'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['bin/load_mols',
             'bin/standardize_mols',
             'bin/peek_hdf',
             'bin/substruct_mols',
             'bin/classify_frags',
             'bin/deglyco_mols',
             'bin/map_frags',
             'bin/chunk_sdf',
             'bin/subset_mols',
             'bin/annotate_pnp',
             'bin/plot_pipeline_results',
             ],
    # package_data={'npfc': ['data/deglyco_mols.knwf']},  # key: package name, value: list of data files
    include_package_data=True,
)
