#!/home/comp/cschaoyang/SOFTWARE/ANACONDA3/V2019/bin/python
import argparse
import os
import sys
import datetime

import pipeline
import utility

def checkFQ(args):
	# check the input fq
	if (not os.path.exists(args.input_fastq1) or not os.path.exists(args.input_fastq1)):
		print("Error: The input fastq 1 ot input fastq 2 is unaccessible. Please check.")
		sys.exit(1)
	elif((args.input_type == "TELLSeq") & (not args.index_fastq)):
		print("Error: The index FQ file is required for TELLSeq sequencing technology.")

	# create output dir
	dirname, filename = os.path.split(os.path.abspath(args.output_fastq1))
	if(os.path.exists(dirname)):
		print("The output directory is:" + dirname)
	else:
		os.makedirs(dirname)

	# start the conversion
	pipeline.moduleFQ(args)


def checkBAM(args):
	# check the input fq
	if (not os.path.exists(args.input_barcoded_fastq1) or not os.path.exists(args.input_barcoded_fastq2)):
		print("Error: The input barcoded fastq 1 ot barcoded fastq 2 is unaccessible. Please check.")
		sys.exit(1)

	# check the reference
	if (not os.path.exists(args.reference)):
		print("Error: The reference genome is unaccessible. Please check.")
		sys.exit(1)
		
	# create output dir
	dirname, filename = os.path.split(os.path.abspath(args.outfile))
	if(os.path.exists(dirname)):
		print("The output directory is: " + dirname)
	else:
		os.makedirs(dirname)

	#tmp
	print(args.read_group)
	# start the alignment
	pipeline.moduleBAM(args)


def checkSNV(args):
	# check the input bam
	if (not os.path.exists(args.bam)):
		print("Error: The input bam file is unaccessible. Please check.")
		sys.exit(1)

	# check the reference
	if (not os.path.exists(args.reference)):
		print("Error: The reference genome is unaccessible. Please check.")
		sys.exit(1)

	# create output dir
	dirname, filename = os.path.split(os.path.abspath(args.outfile))
	if(os.path.exists(dirname)):
		print("The output directory is: " + dirname)
	else:
		os.makedirs(dirname)
	
	# start SNV/INDEL calling
	pipeline.moduleSNV(args)


def checkSV(args):
	# check the input bam
	if (not os.path.exists(args.bam)):
		print("Error: The input bam file is unaccessible. Please check.")
		sys.exit(1)
        
        # check the reference
	if (not os.path.exists(args.reference)):
		print("Error: The reference genome is unaccessible. Please check.")
		sys.exit(1)
	
	# check required database
	if (args.application == "Aquila"):
		if(not os.path.exists(args.uniqness)):
			print("Error: The uniqness database is required to run Aquila. Please check.")
			sys.exit(1)
		if(not os.path.exists(args.vcf)):
			print("Error: The prior SNPs is required to run Aquila. Please check.")
			sys.exit(1)

	elif (args.application == "VALOR"):
		if(not os.path.exists(args.sonic)):
			print("Error: The sonic database is required to run VALOR. Please check.")		
			sys.exit(1)

        # create output directory
	dirname, filename = os.path.split(os.path.abspath(args.outfile))
	if(os.path.exists(dirname)):
		print("The output directory is:" + dirname)
	else:
		os.makedirs(dirname)

        # start SV calling
	pipeline.moduleSV(args)


def checkPHASE(args):
	# check the input vcf
	if (not os.path.exists(args.bam)):
		print("Error: The vcf file is unaccessible. Please check.")
		sys.exit(1)

	# check the input bam
	if (not os.path.exists(args.bam)):
		print("Error: The input bam file is unaccessible. Please check.")
		sys.exit(1)

	# check the reference
	if (not os.path.exists(args.reference)):
		print("Error: The reference genome is unaccessible. Please check.")
		sys.exit(1)
	
	# check the output directory
	dirname, filename = os.path.split(os.path.abspath(args.outfile))
	if(os.path.exists(dirname)):
		print("The output directory is:" + dirname)
	else:
		os.makedirs(dirname)
		
	# start	phasing variants
	pipeline.modulePHASE(args)
	

def checkMKFQ(args):
	# start simulating reads
	pipeline.moduleMKFQ(args)


def check_environment(args):
	# check the dependency
	pipeline.moduleENV(args)

def checkWGS(args):	
	# check the input sample information
	if (not os.path.exists(args.sample_info)):
		print("Error: The sample information file is unaccessible. Please check.")
		sys.exit(1)
	
	# check the database
	# check the reference genome
	dir_ref = args.database + "/" + "GRCH38"
	dir_whitelist = args.database + "/" + "WhiteList"
	dir_uniqness_map = args.database + "/" + "Uniqness_map"
	
	# check the reference genome
	if (not os.path.exists(dir_ref)):
		print("Error: The human genome database is unaccessible. Please check.")
		sys.exit(1)

	# check the whitelist
	if (not os.path.exists(dir_whitelist)):
		print("Error: The whitelist database is unaccessible. Please check.")
		sys.exit(1)

	# check the uniqness_map
	if (not os.path.exists(dir_uniqness_map)):
		print("Error: The uniqness database is unaccessible. Please check.")
		sys.exit(1)

	# create output directory
	if(os.path.exists(args.outdir)):
		print("The output directory is:" + args.outdir)
	else:
		os.makedirs(args.outdir)

	# start the WGS analysis pipeline
	pipeline.moduleWGS(args)


def main():
	dirname, filename = os.path.split(os.path.abspath(sys.argv[0]))
	parser = argparse.ArgumentParser(description="Linked Reads ToolKit",usage="LRTK version 1.0.0\n",conflict_handler='resolve')
	parser.set_defaults(abs=dirname)
	subparsers = parser.add_subparsers()
	
	# checkENV section
	checkENV = subparsers.add_parser("checkENV",
			help="To check the environment")
	checkENV.set_defaults(func=check_environment)

	# MKFQ section
	mkfq = subparsers.add_parser("MKFQ", 
			help="To simulate the linked reads")
	mkfq.add_argument('-CF','--config_file',required=True,
			help='The directory to store config_files')
	mkfq.add_argument('-IT','--input_type',default='stLFR', choices=['10x', 'stLFR', 'TELLSeq'],
			help='The sequencing type for FQs')
	mkfq.set_defaults(func=checkMKFQ)

	#FQCONV
	fqconver = subparsers.add_parser("FQCONVER", help="To convert different FASTQ formats")
	fqconver.add_argument('-I1','--input_fastq1',required=True,
			help='The first file for paired FQs')
	fqconver.add_argument('-I2','--input_fastq2',required=True,
			help='The second file for paired FQs')
	fqconver.add_argument('-ID','--index_fastq',required=False,
			help='The index file for paired FQs')
	fqconver.add_argument('-IT','--input_type',default='stLFR', choices=['10x', 'stLFR', 'TELLSeq'],
			help='The sequencing type for input FQs')
	fqconver.add_argument('-O1','--output_fastq1',required=True,
			help='The first file for output paired FQs')
	fqconver.add_argument('-O2','--output_fastq2',required=True,
			help='The second file for output paired FQs')
	fqconver.add_argument('-OT','--output_type', default='ULRF', choices=['10x', 'ULRF'],
			help='The sequencing type for output FQs')
	fqconver.add_argument('-B', '--barcodes', required=False,
			help='The barcodes list used by different technology')
	fqconver.add_argument('-F','--filter', default='Yes', choices=['Yes', 'No'],
			help='To filt the input FQs')
	fqconver.add_argument('-S','--sort', default='Yes', choices=['Yes', 'No'],
			help='To sort the output FQs')
	fqconver.add_argument('-T','--threads',default=1,
			help='Number of threads')

	fqconver.set_defaults(func=checkFQ)

	#ALIGN
	align = subparsers.add_parser("ALIGN", help="To map reads to reference genome")
	align.add_argument('-BQ1','--input_barcoded_fastq1',required=True,
			help='The first input file for barcoded paired FQs')
	align.add_argument('-BQ2','--input_barcoded_fastq2',required=True,
			help='The second input file for bararcoded paired FQs')
	align.add_argument('-FQ1','--input_fastq1',required=False,
			help='The first input file for no Barcoded paired FQs')
	align.add_argument('-FQ2','--input_fastq2',required=False,
			help='The first input file for no Barcoded paired FQs')
	align.add_argument('-RG','--read_group',default='@RG\tID:example\tSM:example',
			help='The read group string')
	align.add_argument('-R','--reference',required=True,
			help='The reference sequence to align')
	align.add_argument('-O','--outfile',required=True,
			help='The output bamfile')
	align.add_argument('-S','--sort', default="Yes", choices=['Yes', 'No'],
			help='To sort bam file by barcode')
	align.add_argument('-M','--mark_duplication', default='Yes', choices=['Yes', 'No'],
			help='To mark the duplicated reads by barcode')
	align.add_argument('-P','--platform',default='10x', choices=['10x', 'stLFR', 'TELLSeq'],
			help='linked-reads technology')
	align.add_argument('-T','--threads', default=1,
			help='Number of threads')

	align.set_defaults(func=checkBAM)	

	#SNV calling
	snp = subparsers.add_parser("SNV", help="To call small variants such as SNPs and INDELs")
	snp.add_argument('-B','--bam',required=True,
			help='The aligned bam file for variation calling')
	snp.add_argument('-R','--reference',required=True,
			help='The reference sequence')
	snp.add_argument('-A','--application',default='FreeBayes', choices=['FreeBayes', 'Samtools', 'GATK'],
			help='SNV calling tools (FreeBayes, Samtools, GATK)')
	snp.add_argument('-T','--threads',required=False, default=1,
			help='Number of threads')
	snp.add_argument('-O','--outfile',required=True,
			help='output file')
	
	snp.set_defaults(func=checkSNV)
	
	#SV calling
	sv = subparsers.add_parser("SV", help="To call large structural variants")
	sv.add_argument('-B', '--bam', required=True,
			help='The aligned bam file for variation calling')
	sv.add_argument('-V', '--vcf', required=False,
			help='The precalled variants required by Aquila')
	sv.add_argument('-R', '--reference', required=True,
			help='The reference sequence')
	sv.add_argument('-A', '--application', default='Aquila', choices=['Aquila', 'LinkedSV', 'VALOR'],
			help='The available SV calling tools')
	sv.add_argument('-T', '--threads', default=1,
			help='Number of threads')
	sv.add_argument('-U', '--uniqness', required=False,
			help='The uniqness database required by Aquila')
	sv.add_argument('-S', '--sonic', required=False,
			help='The sonic database required by VALOR')
	sv.add_argument('-O', '--outfile', required=True,
			help='The output file')
	
	sv.set_defaults(func=checkSV)

	#PHASE
	phase = subparsers.add_parser("PHASE", help="To phase germline variations")
	phase.add_argument('-B', '--bam', required=True,
			help='The aligned bam file for variation phasing')
	phase.add_argument('-V', '--vcf', required=True,
			help='The detected variants to phase')
	phase.add_argument('-R', '--reference', required=True,
			help='The reference sequence')
	phase.add_argument('-A', '--application', default='HapCUT2', choices=['HapCUT2','WhatsHap','SpecHap'],
			help='The available phasing tools')
	phase.add_argument('-T', '--threads', default=1,
			help='Number of threads')
	phase.add_argument('-O', '--outfile', required=True,
			help='The output file')
	
	phase.set_defaults(func=checkPHASE)

	#WGS
	wgs = subparsers.add_parser("WGS", help="Automatic pipeline to perform SNV/INDEL/SV calling and phasing")
	wgs.add_argument('-SI','--sample_info',required=True,
			help='The descriptive text for input samples (Format:"ID FQ1 FQ2 INDEX Linked_read_technology")')
	wgs.add_argument('-OD','--outdir',required=True,
			help='The output directory')
	wgs.add_argument('-DB','--database',required=True,
			help='The reference sequence to align')
	wgs.add_argument('-RG','--read_group',default='@RG\tID:example\tSM:example',
			help='The read group')
	wgs.add_argument('-T', '--threads', default=1,
			help='Number of threads')
	
	wgs.set_defaults(func=checkWGS)
	
	args = parser.parse_args()
	
	try:
        	args.func(args)
	except AttributeError:
		parser.error("Please try 'LRTK -h' for usage information")

if __name__ == '__main__':
	main()

