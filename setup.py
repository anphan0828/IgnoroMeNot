from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.4'
DESCRIPTION = 'Ignorome finder'

with open("README.md","r") as fh:
    long_description=fh.read()

# Setting up
setup(
        name="ignoromenot",
        version=VERSION,
        author="An Phan",
        author_email="<ahphan@iastate.edu>",
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type = 'text/markdown',
        packages=find_packages(),
        install_requires=[],
        py_modules=["ignoromenot"],
        package_dir={'': 'src'},
        url="https://github.com/anphan0828/IgnoroMeNot.git",
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
        ]
)


