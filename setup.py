#!/usr/bin/env python3
#
# Copyright (C) 2017-2021
#

__author__ = 'Simon Cabello'
__copyright__ = 'Copyright (C) 2021'
__license__ = 'GNU General Public License'
__version__ = '0.2.1'
__email__ = 's-cabelloaguilar@chu-montpellier.fr'
__status__ = 'prod'

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ifCNV",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="ifCNV",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SimCab-CHU/ifCNV",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    install_requires=[
        'numpy>=1.21',
        'pandas>=1.3',
        'scikit-learn>=1.0.1',
        'plotly>=5.4',
        "pybedtools>=0.8.2"
    ],
    scripts=['script/ifCNV'],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/SimCab-CHU/ifCNV/issues',
        'Source': 'https://github.com/SimCab-CHU/ifCNV',
    }
)
