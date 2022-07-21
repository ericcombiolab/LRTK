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
    version='1.0',
    description='Linked-Read ToolKit',
    long_description='Linked-Read ToolKit',
    long_description_content_type='text/markdown',
    author='YANG CHAO',
    author_email='yangchaogab@gmail.com',
    url='https://github.com/ericcombiolab/LRTK',
    packages=['lrtk'],
    install_requires=['pysam','pandas','bwa','picard','samtools','freebayes','hapcut2', 'aquila', 'vcflib'],
    entry_points={'console_scripts':['entry=packagename.scriptname:function',],
    include_package_data=True,
    license='MIT',
    zip_safe=False
)
