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
    version='1.9',
    description='Linked-Read ToolKit',
    author='YANG CHAO',
    author_email='yangchaogab@gmail.com',
    packages=find_packages(),
    install_requires=['numpy','pysam', "pytorch", 'scipy', 'scikit-learn','sortedcontainers', 'aquila', 'bam-readcount','bcftools', 'bedtools','bwa','fastp', 'flye', 'freebayes','gatk','hapcut2','jellyfish', 'megahit', 'metabat2', 'parallel','picard', 'pigz', 'quickmerge', 'r-base', 'r-ggplot2', 'r-ggforce', 'r-vegan', 'r-factoextra', 'r-proxy', 'r-clusterSim', 'samtools', 'snakemake', 'seqtk', 'spades', 'whatshap','vcflib'],
    package_data={'script' : ['*pl',
                            '*R',
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
                            'Pangaea/*',
                            'Pangaea/bin/*',
                            'Pangaea/bin/Lathe/*',
                            'Pangaea/bin/Lathe/scripts/*',
                            'Pangaea/cpptools/*',
                            'Pangaea/cpptools/lib/*',
                            'Pangaea/cpptools/include/*',
                            'Pangaea/cpptools/include/ThreadPool/*',
                            'Pangaea/cpptools/include/htslib/*',
                            'Pangaea/models/*',
                             ]},
    entry_points={'console_scripts':['lrtk=script.lrtk:main']},
    license='MIT',
    zip_safe=False
)
