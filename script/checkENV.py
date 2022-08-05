#!/usr/bin/python

import argparse
import os
import gzip
import sys
import shutil
from utility import *

def main():
	##Here we will check the requirements about dependency
	##check SV calling tools: Aquila
	tool_path = shutil.which("Aquila_step1")
	if tool_path is None:
		print('Error:Aquila is not installed or included in the environment ...')
		print('using conda to install aquila ...')
		cmd="conda install -c bioconda aquila"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing aquila ...')


	##check the alignment tool: bwa
	tool_path = shutil.which("bwa")
	if tool_path is None:
		print('Error: bwa is not installed or included in the environment ...')
		print('using conda to install bwa ...')
		cmd="conda install -c bioconda bwa"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing bwa ...')
	
	##check the SNV calling tool:freebayes
	tool_path = shutil.which("freebayes")
	if tool_path is None:
		print('Error: freebayes is not installed or included in the environment ...')
		print('using conda to install freebayes ...')
		cmd="conda install -c bioconda freebayes"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing freebayes ...')

	##check the SNV calling tool:GATK
	tool_path = shutil.which("gatk3")
	if tool_path is None:
		print('Error: GATK is not installed or included in the environment ...')
		print('using conda to install gatk ...')
		cmd="conda install -c bioconda gatk"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing GTAK ...')

	##check the phasing tool:HapCUT2
	tool_path = shutil.which("hapcut2")
	if tool_path is None:
		print('Error: HapCUT2 is not installed or included in the environment ...')
		print('using conda to install hapcut2 ...')
		cmd="conda install -c bioconda hapcut2"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing HapCUT2 ...')

	##check the SV calling tool:LinkedSV
	tool_path = shutil.which("linkedsv.py")
	if tool_path is None:
		print('Error: LinkedSV is not installed or included in the environment ...')
		print('Please install LinkedSV from source ...')
#		cmd="git clone https://github.com/WGLab/LinkedSV.git "
#		subprocess.call(cmd, shell=True)
#		cmd="cd LinkedSV && sh build.sh"
#		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing LinkedSV ...')

	##check the tool: parallel
	tool_path = shutil.which("parallel")
	if tool_path is None:
		print('Error: parallel is not installed or included in the environment ...')
		print('using conda to install parallel ...')
		cmd="conda install -c bioconda parallel"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing parallel ...')

	##check the bam processing tool: picard
	tool_path = shutil.which("picard")
	if tool_path is None:
		print('Error: picard is not installed or included in the environment ...')
		print('using conda to install picard ...')
		cmd="conda install -c bioconda picard"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing picard ...')
	
	##check the bam processing tool:samtools
	tool_path = shutil.which("samtools")
	if tool_path is None:
		print('Error: samtools is not installed or included in the environment ...')
		print('using conda to install samtools ...')
		cmd="conda install -c bioconda samtools"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing samtools ...')

	##check the phasing tool:SpecHap
	tool_path = shutil.which("SpecHap")
	if tool_path is None:
		print('Error: SpecHap is not installed or included in the environment ...')
		print('Please install SpecHap from source ...')
#		cmd="git clone https://github.com/deepomicslab/SpecHap.git"
#		subprocess.call(cmd, shell=True)
#		cmd="cd SpecHap && mkdir build && cd build && cmake .. && make && make install"
#		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing SpecHap ...')

	##check the phasing tool: whatshap
	tool_path = shutil.which("whatshap")
	if tool_path is None:
		print('Error: whatshap is not installed or included in the environment ...')
		print('using conda to install whatshap ...')
		cmd="conda install -c bioconda whatshap"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing whatshap ...')
	
	##check the SV calling tool: valor
	tool_path = shutil.which("valor")
	if tool_path is None:
		print('Error: valor is not installed or included in the environment ...')
		print('Please install VALOR from source ...')
#		cmd="git clone https://github.com/BilkentCompGen/valor.git --recursive"
#		subprocess.call(cmd, shell=True)
#		cmd="cd valor && make libs && make && cp -r ../valor ymp"
#		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing valor ...')

	tool_path = shutil.which("vcfsort")
	if tool_path is None:
		print('Error: valor is not installed or included in the environment ...')
		print('using conda to install vcflib ...')
		cmd="conda install -c bioconda vcflib"
		subprocess.call(cmd, shell=True)
	else:
		print('Pass: using existing vcflib ...')

if __name__ == "__main__":
	main()
