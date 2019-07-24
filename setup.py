import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


# begin

setuptools.setup(
    name="npfc",
    version='0.3.8-1-g4b59d9a',
    author="Stephane Bourg, Jose-Manuel Gally",
    author_email="stephane.bourg@crns.fr",
    description="A package for classifying fragment combinations in molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=[
          'pip',
      ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    keywords=['chemoinformatics', 'natural products', 'fragments', 'chemical biology'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=[
             'bin/chunk_sdf',
             'bin/load_mols',
             'bin/deglyco_mols',
             'bin/standardize_mols',
             'bin/filter_dupl_mols',
             'bin/murcko_mols',
             'bin/subset_mols',
             'bin/substruct_mols',
             'bin/classify_frags',
             'bin/map_frags',
             'bin/annotate_pnp',
             'bin/compute2D_mols',
             'bin/plot_pipeline_results',
             'bin/peek_hdf',
             ],
    # package_data={'npfc': ['data/deglyco_mols.knwf']},  # key: package name, value: list of data files
    include_package_data=True,
)
