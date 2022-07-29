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
    packages=find_packages(),
    install_requires=['pysam','pandas','bwa','picard','samtools','freebayes','hapcut2', 'aquila', 'vcflib'],
    package_data={'lrtk' : ['database/*',
    		            'example/*',
                            'example/FQs/*',
                            'example/FQs/simulation/resource/*',
                            'example/FQs/simulation/diploid_config*',
                            'example/LargeFQs/*',
                            'script/*pl*,
                            'script/correct_barcode*',
                            'script/EMA/*',
                            'script/EMA/bwa/*',
                            'script/EMA/cpp/*',
                            'script/EMA/include/*',
                            'script/EMA/obj/*',
                            'script/EMA/src/*',
                            'script/header/*',
                            'script/long_fragment/*',
                            'script/simulation/resource/*',
                            'script/simulation/diploid_config/*',
                            'script/LinkedSV/*',
                            'script/valor/*',
                            'script/SpecHap/*'
                             ]},
    entry_points={'console_scripts':['lrtk=lrtk.script.lrtk:main']},
    license='MIT',
    zip_safe=False
)
