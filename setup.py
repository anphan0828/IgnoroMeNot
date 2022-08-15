from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.2'

# Setting up
setup(
        name="ignoromenot",
        version=VERSION,
        author="An Phan",
        author_email="<ahphan@iastate.edu>",
        description="Finding ignorome",
        long_description="This program outputs a list of ignorome genes highly associated with other well-annotated genes",
        readme="README.md",
        packages=find_packages(),
        install_requires=[],
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
        ]
)


