#!/home/comp/cschaoyang/SOFTWARE/ANACONDA3/V2019/bin/python
import argparse
import os
import sys
import datetime

sys.path.append(os.path.dirname(__file__))

from pipeline import *
from utility import *

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
		print("The output directory is: " + dirname)
	else:
		os.makedirs(dirname)

	# start the conversion
	moduleFQ(args)


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
	moduleBAM(args)


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
	moduleSNV(args)


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
	moduleSV(args)


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
	modulePHASE(args)
	

def checkMKFQ(args):
	# start simulating reads
	moduleMKFQ(args)


def check_environment(args):
	# check the dependency
	moduleENV(args)

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
	moduleWGS(args)


def main():
	parser = argparse.ArgumentParser(description="Linked Reads ToolKit",usage="lrtk version 1.4\n",conflict_handler='resolve')
	subparsers = parser.add_subparsers()
	
	# checkENV section
	checkENV = subparsers.add_parser("checkENV",
			help="Check the working environment")
	checkENV.set_defaults(func=check_environment)

	# MKFQ section
	mkfq = subparsers.add_parser("MKFQ", 
			help="Simulate linked-reads")
	mkfq.add_argument('-CF','--config_file',required=True,
			help='The path to config_files for simulation')
	mkfq.add_argument('-IT','--input_type',default='stLFR', choices=['10x', 'stLFR'],
			help='Input sequencing technology. Users can choose from (10x,stLFR)')
	mkfq.set_defaults(func=checkMKFQ)

	#FQCONV
	fqconver = subparsers.add_parser("FQCONVER", help="Convert FASTQ formats")
	fqconver.add_argument('-I1','--input_fastq1',required=True,
			help='Input fastq file (uncompressed FASTQ format), the first read of paired linked-read sequencing data')
	fqconver.add_argument('-I2','--input_fastq2',required=True,
			help='Input fastq file (uncompressed FASTQ format), the second read of paired linked-read sequencing data')
	fqconver.add_argument('-ID','--index_fastq',required=False,
			help='Input index file (uncompressed FASTQ format) for paired linked-read sequencing data.')
	fqconver.add_argument('-IT','--input_type',default='stLFR', choices=['10x', 'stLFR', 'TELLSeq'],
			help='Input sequencing technology. Users can choose from (10x,stLFR,TELLSeq).')
	fqconver.add_argument('-O1','--output_fastq1',required=True,
			help='Output fastq file, the first read of paired linked-read sequencing data')
	fqconver.add_argument('-O2','--output_fastq2',required=True,
			help='Output fastq file, the second read of paired linked-read sequencing data.')
	fqconver.add_argument('-OT','--output_type', default='ULRF', choices=['ULRF'],
			help='Output fastq format. The unified linked read format (ULRF) is recommended')
	fqconver.add_argument('-B', '--barcodes', required=False,
			help='The reference barcode whitelist files for 10x and stLFR technologies')
	fqconver.add_argument('-F','--filter', default='Yes', choices=['Yes', 'No'],
			help='Users can choose from (Yes, No). "Yes" indicates that LRTK will use fastp to filter reads.')
	fqconver.add_argument('-S','--sort', default='Yes', choices=['Yes', 'No'],
			help='Users can choose from (Yes, No). "Yes" indicates that LRTK will sort the reads based on barcodes.')
	fqconver.add_argument('-T','--threads',default=1,
			help='Number of threads, this determines the number of threads used for bwa and samtools')

	fqconver.set_defaults(func=checkFQ)

	#ALIGN
	align = subparsers.add_parser("ALIGN", help="Align reads to the reference genome")
	align.add_argument('-BQ1','--input_barcoded_fastq1',required=True,
			help='Input fastq file (uncompressed FASTQ format), the first read of paired linked-read sequencing data (with barcode).')
	align.add_argument('-BQ2','--input_barcoded_fastq2',required=True,
			help='Input fastq file (uncompressed FASTQ format), the second read of paired linked-read sequencing data (with barcode).')
	align.add_argument('-FQ1','--input_fastq1',required=False,
			help='Input fastq file (uncompressed FASTQ format), the first read of paired linked-read sequencing data (without barcode).')
	align.add_argument('-FQ2','--input_fastq2',required=False,
			help='Input fastq file (uncompressed FASTQ format), the second read of paired linked-read sequencing data (without barcode).')
	align.add_argument('-RG','--read_group',default='@RG\tID:example\tSM:example',
			help='Full read group string (e.g. @RG\tID:foo\tSM:bar)')
	align.add_argument('-R','--reference',required=True,
			help='The indexed human reference genome file')
	align.add_argument('-O','--outfile',required=True,
			help='The output alignment file.')
	align.add_argument('-S','--sort', default="Yes", choices=['Yes', 'No'],
			help='Users can choose from (Yes, No). "Yes" means LRTK will use samtools to sort alignment files based on genomic coordinate.')
	align.add_argument('-M','--mark_duplication', default='Yes', choices=['Yes', 'No'],
			help='Users can choose from (Yes, No). "Yes" means LRTK will use picard to mark the duplicated reads using barcode information')
	align.add_argument('-P','--platform',default='10x', choices=['10x', 'stLFR', 'TELLSeq'],
			help='Input sequencing technology. Users can choose from (10x,stLFR,TELLSeq).')
	align.add_argument('-T','--threads', default=1,
			help='Number of threads, this determines the number of threads used for ema, bwa and samtools.')

	align.set_defaults(func=checkBAM)	

	#SNV calling
	snp = subparsers.add_parser("SNV", help="Detect SNVs and INDELs")
	snp.add_argument('-B','--bam',required=True,
			help='The alignment file (.bam).')
	snp.add_argument('-R','--reference',required=True,
			help='The indexed human reference genome file.')
	snp.add_argument('-A','--application',default='FreeBayes', choices=['FreeBayes', 'Samtools', 'GATK'],
			help='The SNV/INDEL caller. Users can choose from (FreeBayes, Samtools, GATK).')
	snp.add_argument('-T','--threads',required=False, default=1,
			help='Number of threads, this determines the number of threads used for SNV/INDEL caller.')
	snp.add_argument('-O','--outfile',required=True,
			help='The final VCF file to write.')
	
	snp.set_defaults(func=checkSNV)
	
	#SV calling
	sv = subparsers.add_parser("SV", help="Detect structural variations")
	sv.add_argument('-B', '--bam', required=True,
			help='The alignment file (.bam).')
	sv.add_argument('-V', '--vcf', required=False,
			help='The precalled SNV/INDEL variants required for Aquila')
	sv.add_argument('-R', '--reference', required=True,
			help='The indexed human reference genome file.')
	sv.add_argument('-A', '--application', default='Aquila', choices=['Aquila', 'LinkedSV', 'VALOR'],
			help='The SV caller. Users can choose from (Aquila, LinkedSV, VALOR).')
	sv.add_argument('-T', '--threads', default=1,
			help='Number of threads,this determines the number of threads used for SV caller.')
	sv.add_argument('-U', '--uniqness', required=False,
			help='The uniqness database is required for Aquila.')
	sv.add_argument('-S', '--sonic', required=False,
			help='The sonic database is required for VALOR.')
	sv.add_argument('-O', '--outfile', required=True,
			help='The final VCF file to write.')
	
	sv.set_defaults(func=checkSV)

	#PHASE
	phase = subparsers.add_parser("PHASE", help="Phase germline variations")
	phase.add_argument('-B', '--bam', required=True,
			help=' The alignment file (.bam).')
	phase.add_argument('-V', '--vcf', required=True,
			help='The detected variants to phase')
	phase.add_argument('-R', '--reference', required=True,
			help='The indexed human reference genome file.')
	phase.add_argument('-A', '--application', default='HapCUT2', choices=['HapCUT2','WhatsHap','SpecHap'],
			help='The variant phasing tool. Users can choose from (HapCUT2, WhatsHap, SpecHap).')
	phase.add_argument('-T', '--threads', default=1,
			help='Number of threads, this determines the number of threads used for phasing tools.')
	phase.add_argument('-O', '--outfile', required=True,
			help='The final phased VCF file to write.')
	
	phase.set_defaults(func=checkPHASE)

	#WGS
	wgs = subparsers.add_parser("WGS", help="Run the whole pipeline")
	wgs.add_argument('-SI','--sample_info',required=True,
			help='The path to input sample information file. (sample_info Format:"ID FQ1 FQ2 INDEX Linked_read_technology")')
	wgs.add_argument('-OD','--outdir',required=True,
			help='The output directory')
	wgs.add_argument('-DB','--database',required=True,
			help='The default database containing reference genome and barcode whitelist file')
	wgs.add_argument('-RG','--read_group',default='@RG\tID:example\tSM:example',
			help='Full read group string (e.g. @RG\tID:foo\tSM:bar).')
	wgs.add_argument('-T', '--threads', default=1,
			help='Number of threads, this determines the number of threads used for alignment and variats detection and phasing.')
	
	wgs.set_defaults(func=checkWGS)
	
	args = parser.parse_args()
	
	try:
        	args.func(args)
	except AttributeError:
		parser.error("Please try 'lrtk -h' for usage information")

if __name__ == '__main__':
	main()

