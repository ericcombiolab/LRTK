import os
import sys
import subprocess

script_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(script_path)

from . import conversion
from . import alignment
from . import small_variants
from . import large_variants
from . import phasing
#from . import 
from utility import *


def moduleFQ(args):	
	# To convert different fastq formats
	# convert 10x to ULRF
	if((args.input_type == "10x") & (args.output_type == "ULRF")): #convert 10x to ULRF
		conversion.TENx2ULRF(args.input_fastq1, args.input_fastq2, args.output_fastq1, args.output_fastq2, args.barcodes, args.filter, args.sort, args.threads, args.abs)
	elif((args.input_type == "stLFR") & (args.output_type == "ULRF")): # convert stLFR to ULRF
		conversion.stLFR2ULRF(args.input_fastq1, args.input_fastq2, args.output_fastq1, args.output_fastq2, args.barcodes, args.filter, args.sort, args.threads, args.abs)
	elif((args.input_type == "TELLSeq") & (args.output_type == "ULRF")): # convert TELL-Seq to ULRF
		conversion.TELLSeq2ULRF(args.input_fastq1, args.input_fastq2, args.index_fastq, args.output_fastq1, args.output_fastq2, args.filter, args.sort, args.threads, args.abs)
	elif((args.input_type == "stLFR") & (args.output_type == "10x")): # convert stLFR to 10x
		conversion.stLFR210x(args.input_fastq1, args.input_fastq2, args.output_fastq1, args.output_fastq2, args.abs)
	elif((args.input_type == "TELLSeq") & (args.output_type == "10x")): # convert TELL-Seq to 10x 
		conversion.TELLSeq210x(args.input_fastq1, args.input_fastq2, args.index_fastq, args.output_fastq1, args.output_fastq2, args.barcodes, args.abs)

def moduleMKFQ(args):
	# To simulate linked reads
	# for 10x
	if(args.input_type == "10x"):
		run_cmd([
			"python", 
			script_path+"simulate_reads.py",
			args.config_file,
			],
			"simulate_read"
			)
	# for stlfr
	elif(args.input_type == "stLFR"):
		run_cmd([
			"python",
			script_path+"simulate_reads_stLFR.py",
			args.config_file,
			],
			"simulate_read"
			)

def moduleBAM(args):
	if((args.platform == "10x")):
		alignment.align_10x(args.input_barcoded_fastq1, args.input_barcoded_fastq2, args.input_fastq1, args.input_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.threads, args.abs)
	elif(args.platform == "stLFR"):
		alignment.align_stLFR(args.input_barcoded_fastq1, args.input_barcoded_fastq2, args.input_fastq1, args.input_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.threads, args.abs)
	elif(args.platform == "TELLSeq"):
		alignment.align_TELLSeq(args.input_barcoded_fastq1, args.input_barcoded_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.threads, args.abs)
	else:
		print("Error: unknown linked-read technology!")


def moduleSNV(args):	
	if((args.application == "FreeBayes")):
		SmallVariants.calling_freebayes(args.bam, args.reference, args.threads, args.outfile, args.abs)
	elif(args.platform == "Samtools"):
		SmallVariants.calling_samtools(args.bam, args.reference, args.threads, args.outfile, args.abs)
		print("SNV calling with Samtools\n")
	elif(args.platform == "GATK"):
		print("SNV calling with GATK\n")
	else:
		print("Error: unknown SNV/INDEL caller!")


def moduleSV(args):
	if(args.application == "Aquila"):
		print("Here is moduleSV")
		LargeVariants.calling_Aquila(args.bam, args.vcf, args.reference, args.uniqness, args.threads, args.outfile, args.abs)
	elif(args.application == "LinkedSV"):
		LargeVariants.calling_LinkedSV(args.bam, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "VALOR"):
		LargeVariants.calling_LinkedSV(args.bam, args.reference, args.sonic, args.threads, args.outfile, args.abs)
	else:
		print("Error: unknown SV caller!")


def modulePHASE(args):
	###make preparations
	if((args.application == "HapCUT2")):
		haplotyping.Phasing_HAPCUT2(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "WhatsHap"):
		haplotyping.Phasing_WhatsHap(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "SpecHap"):
		haplotyping.Phasing_SpecHap(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)
	else:
		print("Error: unknown phasing tool!")


def moduleWGS(args):
	###whole pipeline
	###FASTQ processing && mapping
	sample_infile = args.sample_info
	INFILE=open(sample_infile,"rt")
	line=INFILE.readline()
	while(line):
		if not line.startwith("#"):
			arr=line.strip().split("\t")
			barcode = arr[0]
			input_fastq1 = arr[1]
			input_fastq2 = arr[2]
			index_fastq  = arr[3]
			tech         = arr[4]

			whitelist=args.index_fastq
			threads=args.threads
			outdir=args.outdir
			ref=args.reference

			outdirfq  = outdir + "/FASTQ"
			outdirbam = outdir + "/BAM"
			outfq1 = outdirfq + "Clean.R1.fq"
			outfq2 = outdirfq + "Clean.R2.fq"
			wbfq1  = outdirfq + "Clean.R1.fq"
			wbfq2  = outdirfq + "Clean.R2.fq"
			wobfq1 = outdirfq + "Clean.R1.fq"
			wobfq2 = outdirfq + "Clean.R2.fq"
			bamfile = outdirbam+"Clean.sort.rmdup.bam"
			
			if((tech == "10x")):
				conversion.TENx2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, "Yes", "Yes", threads, bin)
				alignment.align_10x(wbfq1, wbfq2, wobfq1, wobfq2, args.read_group, ref, bamfile, "Yes", "Yes", threads, bin)
				
			elif(tech == "stLFR"):
				conversion.stLFR2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, "Yes", "Yes", threads, bin)
				alignment.align_stLFR(input_fastq1, input_fastq2, outfq1, outfq2, args.read_group, ref, bamfile, "Yes", "Yes", threads, bin)
			elif(tech == "TELLSeq"):
				conversion.TELLSeq2ULRF(input_fastq1, input_fastq2, index_fastq, outfq1, outfq2, "Yes", "Yes", threads, bin)
				alignment.align_TELLSeq(wbfq1, wbfq2, args.read_group, ref, bamfile, "Yes", "Yes", threads, bin)

			else:
				print("Unknown linked-read technology!")

			###Variants calling && phasing
			SmallVariants.calling_freebayes(args.bam, args.reference, args.threads, args.outfile, args.abs)
			LargeVariants.calling_Aquila(args.bam, args.vcf, args.reference, args.uniqness, args.threads, args.outfile, args.abs)
			haplotyping.Phasing_HAPCUT2(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)

		line=INFILE.readline()

	
