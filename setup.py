#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
import sys
from shutil import rmtree
from setuptools import find_packages, setup, Command

# Where the magic happens:
from setuptools import setup, find_packages, Extension, Command
setup(
    name='lrtk',
    version='1.5',
    description='Linked-Read ToolKit',
    author='YANG CHAO',
    author_email='yangchaogab@gmail.com',
    packages=find_packages(),
    install_requires=['numpy','pysam','scipy','sortedcontainers','aquila','bcftools','bwa','fastp','freebayes','gatk','hapcut2','parallel','picard','samtools','whatshap','vcflib'],
    package_data={'script' : ['*pl',
                            '*cpp',
                            '*stlfr',
                            'EMA/*',
                            'EMA/bwa/*',
                            'EMA/cpp/*',
                            'EMA/include/*',
                            'EMA/obj/*',
                            'EMA/src/*',
                            'header/*',
                            'long_fragment/*',
                            'LinkedSV/*',
                            'valor/*',
                            'SpecHap/*',
                             ]},
    entry_points={'console_scripts':['lrtk=script.lrtk:main']},
    license='MIT',
    zip_safe=False
)
