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
    version='1.1',
    description='Linked-Read ToolKit',
    author='YANG CHAO',
    author_email='yangchaogab@gmail.com',
    packages=['lrtk'],
    install_requires=['pysam','pandas','bwa','picard','samtools','freebayes','hapcut2', 'aquila', 'vcflib'],
    entry_points={'console_scripts':['lrtk=lrtk.lrtk:main']},
    license='MIT',
    zip_safe=False
)
