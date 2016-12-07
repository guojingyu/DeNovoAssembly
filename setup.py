import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "de_novo_assembly",
    version = "0.0.1",
    author = "Jingyu Guo",
    description = ("A package assembles short DNA reads into contigs/gaps."),
    license = "MIT",
    keywords = "example documentation tutorial",
    url = "http://packages.python.org/de_novo_assembly",
    packages=['de_novo_assembly', 'tests', 'docs', 'data'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)
