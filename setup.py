import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biobb_dyna",
    version="5.0.0",
    author="Biobb developers",
    author_email="pieter.zanders@bsc.es",
    description="biobb_dyna is the Biobb module collection to create and train ML & DL models.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="Bioinformatics Workflows BioExcel Compatibility",
    url="https://github.com/bioexcel/biobb_dyna",
    project_urls={
        "Documentation": "http://biobb-autoencoders.readthedocs.io/en/latest/",
        "Bioexcel": "https://bioexcel.eu/",
    },
    packages=setuptools.find_packages(exclude=["docs", "test"]),
    package_data={"biobb_dyna": ["py.typed", "json_schemas/*.json"]},
    install_requires=["biobb_common==5.0.0", "mdigest", "networkx", "MDAnalysis"],
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "dccm = biobb_dyna.dyna.dyncorr.dccm:main",
            "create_graph = biobb_dyna.dyna.network.create_graph:main",
            "communities = biobb_dyna.dyna.analysis.communities:main",
            "shortest_paths = biobb_dyna.dyna.analysis.shortest_paths:main",
        ]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
    ],
)
