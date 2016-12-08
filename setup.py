
from setuptools import setup, find_packages

setup(
    name="DeNovoAssembly",
    version="0.0.1",
    packages=find_packages(),
    scripts=['run.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['python == 2.7.11',
                      'nose == 1.3.7',
                      'numpy == 1.11.1',
                      'biopython == 1.68',
                      'networkx == 1.11',
                      'matplotlib == 1.5.3'],

    package_data={
        # If any package contains *.md or *.text, and *.png files, include them:
        '': ['*.txt', '*.md'],
        # as well LICENSE
        '': ['LICENSE'],
        # And include any *.sh files found in the 'example' package, too:
        'example': ['*.sh'],
        # And include any *.fasta files found in the 'data' package, too:
        'data': ['*.fasta'],
    },

    # metadata for upload to PyPI
    author="Jingyu Guo",
    description="This is a package containing an implementation of de novo "
                "assembly algorithm based on de Bruijn graph.",
    license="MIT",
    keywords="de novo assembly, de Bruijn graph, Eulerian",
    url="https://github.com/guojingyu/DeNovoAssembly",
)
