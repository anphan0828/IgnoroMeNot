from setuptools import setup, find_packages
import codecs
import os
from os import path

here = path.abspath(path.dirname(__file__))

VERSION = '0.0.9'
DESCRIPTION = 'Ignorome finder'

with open(path.join(here,"README.md"),encoding='utf-8') as fh:
    long_description=fh.read()

# Setting up
setup(
        name='ignoromenot',
        version=VERSION,
        author="An Phan",
        author_email="<ahphan@iastate.edu>",
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type = 'text/markdown',
        license='GPLv3',
        packages=find_packages(),
        install_requires=['networkx', 'matplotlib','numpy','Biopython','requests','pandas'],
        include_package_data=True,
        # package_dir={'': 'src'},
        entry_points={
            'console_scripts':[
                'ignoromenot = ignoromenot.ignoromenot:main',
            ]
        },
        url="https://github.com/anphan0828/IgnoroMeNot.git",
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
        ]
)


