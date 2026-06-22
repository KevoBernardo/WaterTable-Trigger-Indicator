# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

# The directory containing this file
HERE = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open(path.join(HERE, 'requirements.txt')) as f:
    requirements = f.read().splitlines()

# This call to setup() does all the work
setup(
    name="WaterTableTrigger",
    version="0.0.1",
    description="A simple Python tool to calculate an indicator for water table raise triggers.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    author="Kevin Bernardo",
    author_email="",
    license="SRK/LAME copyright",
    classifiers=[
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent"
    ],
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    package_data={"": ["*.toml"]},
    include_package_data=True,
    install_requires=requirements
)
