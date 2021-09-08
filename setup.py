#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Dr. Pavel Sidorov <pavel.o.sidorov@gmail.com>
#  Copyright 2021 Dr. Timur Gimadiev <timur.gimadiev@gmail.com>
#  This file is part of HPLCAnalysis project.
#
#  HPLCAnalysis is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from pathlib import Path
from setuptools import setup, find_packages


version = '0.0.1'

setup(
    name='HPLCAnalysis',
    version=version,
    packages=find_packages(),
    url='https://github.com/POSidorov/HPLCAnalysis',
    license='LGPLv3',
    author=['Dr. Pavel Sidorov'],
    author_email=['pavel.o.sidorov@gmail.com'],
    python_requires='>=3.9.0',
    install_requires=['pandas>=1.3', 'numpy>=1.21', 'scipy>=1.7', 'plotly>=5.3', 'lmfit>=1.0'],
    long_description=(Path(__file__).parent / 'README.md').open().read(),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Software Development :: Libraries :: Python Modules',
                 'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.9',
                 ]
)
