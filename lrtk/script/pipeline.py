import os
import shutil
import sys
import subprocess

script_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(script_path)

import conversion
import alignment
import small_variants
import large_variants
import phasing 
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


def moduleENV(args):
	#To check the working environment
	cmd="python " + script_path+"checkENV.py"
	subprocess.call(cmd, shell=True)


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
		small_variants.calling_freebayes(args.bam, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "Samtools"):
		small_variants.calling_samtools(args.bam, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "GATK"):
		small_variants.calling_GATK(args.bam, args.reference, args.threads, args.outfile, args.abs)
	else:
		print("Error: unknown SNV/INDEL caller!")


def moduleSV(args):
	if(args.application == "Aquila"):
		print("Here is moduleSV")
		large_variants.calling_Aquila(args.bam, args.vcf, args.reference, args.uniqness, args.threads, args.outfile, args.abs)
	elif(args.application == "LinkedSV"):
		large_variants.calling_LinkedSV(args.bam, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "VALOR"):
		large_variants.calling_VALOR(args.bam, args.reference, args.sonic, args.threads, args.outfile, args.abs)
	else:
		print("Error: unknown SV caller!")


def modulePHASE(args):
	###make preparations
	if((args.application == "HapCUT2")):
		phasing.Phasing_HAPCUT2(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "WhatsHap"):
		phasing.Phasing_WhatsHap(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)
	elif(args.application == "SpecHap"):
		phasing.Phasing_SpecHap(args.bam, args.vcf, args.reference, args.threads, args.outfile, args.abs)
	else:
		print("Error: unknown phasing tool!")


def moduleWGS(args):
	###FASTQ processing && mapping
	sample_infile = args.sample_info
	rg = args.read_group
	INFILE = open(sample_infile,"rt")
	line = INFILE.readline()
	bin = script_path + "/.."
	while(line):
		if line.startswith("#"):
			print("")
		else:
			arr=line.strip().split("\t")
			
			sample_name  = arr[0]
			input_fastq1 = arr[1]
			input_fastq2 = arr[2]
			index_fastq  = arr[3]
			tech         = arr[4]

			threads=args.threads
			ref=args.database + "/GRCH38/genome.fa"
			uniqness=args.database + "/Uniqness_map/"
			outdir=args.outdir
			if args.outdir.startswith("/"):
				outdir = args.outdir
			else:
				outdir = os.getcwd() + "/" + args.outdir
			print(outdir)
			
			outdir_sample = outdir + "/" + sample_name
			resFQ=outdir_sample+"/"+"FQ"
			resBAM=outdir_sample+"/"+"BAM"
			resQC=outdir_sample+"/"+"QC"
			resSNV=outdir_sample+"/"+"SNV"
			resSV=outdir_sample+"/"+"SV"
			resPS=outdir_sample+"/"+"Phasing"
			resRP=outdir_sample+"/"+"Report"

			os.makedirs(outdir_sample)
			os.makedirs(resFQ)
			os.makedirs(resBAM)
			os.makedirs(resQC)
			os.makedirs(resSNV)
			os.makedirs(resSV)
			os.makedirs(resPS)
			os.makedirs(resRP)

			outdirfq  = resFQ + "/" + sample_name
			outdirbam = resBAM + "/" + sample_name
			outfq1 = outdirfq + ".Clean.R1.fq"
			outfq2 = outdirfq + ".Clean.R2.fq"
			wbfq1  = outdirfq + ".Clean.R1.fq.sort.wb.fq"
			wbfq2  = outdirfq + ".Clean.R2.fq.sort.wb.fq"
			wobfq1 = outdirfq + ".Clean.R1.fq.sort.wob.fq"
			wobfq2 = outdirfq + ".Clean.R2.fq.sort.wob.fq"
			QCFQ   = resFQ + "/" + "FASTQ.QC.json"
			
			bamfile = outdirbam + ".Clean.sort.rmdup.bam"
			bamstat = outdirbam + ".Clean.sort.rmdup.bam.AlignmentStat.xls"
			QCBAM   = resQC + "/" + "Fragment.json"
		
			outvcf1 = resSNV + "/" + sample_name+".smallvariants.vcf"
			outvcf2 = resSV  + "/" + sample_name+".largevariants.vcf"
			outvcf3 = resPS  + "/" + sample_name+".phased.vcf"
			QCVCF   = resQC  + "/" + "VCF.json"

			htmlfile = resRP + "/" + sample_name + ".report.html"

			if((tech == "10x")):
				whitelist=args.database + "/WhiteList/white_list_10x_barcode.fa"
				conversion.TENx2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, "Yes", "Yes", threads, bin)
				alignment.align_10x(wbfq1, wbfq2, wobfq1, wobfq2, rg, ref, bamfile, "Yes", "Yes", threads, script_path)
				
			elif(tech == "stLFR"):
				whitelist=args.database + "/WhiteList/white_list_stlfr_barcode.fa"
				conversion.stLFR2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, "Yes", "Yes", threads, bin)
				alignment.align_stLFR(wbfq1, wbfq2, wobfq1, wobfq2, rg, ref, bamfile, "Yes", "Yes", threads, script_path)

			elif(tech == "TELLSeq"):
				
				conversion.TELLSeq2ULRF(input_fastq1, input_fastq2, index_fastq, outfq1, outfq2, "Yes", "Yes", threads, bin)
				alignment.align_TELLSeq(wbfq1, wbfq2, rg, ref, bamfile, "Yes", "Yes", threads, script_path)

			else:
				print("Unknown linked-read technology!")

			samtools_path = shutil.which("samtools")
			cmd = "python " + script_path+ "/Fragment2json.v2.py" + " -b " + bamfile + " -f " + resBAM + "/constructedfragment.xls" + " -i " + bamfile + ".insert.list" + " -s " + resBAM + "/constructedfragment.stat.xls" + " -t " + samtools_path + " -o " + QCBAM
			print(cmd)
			subprocess.call(cmd, shell=True)

			###Variants calling && phasing
			small_variants.calling_freebayes(bamfile, ref, threads, outvcf1, bin)
			phasing.Phasing_HAPCUT2(bamfile, outvcf1, ref, threads, outvcf3, bin)
			large_variants.calling_Aquila(bamfile, outvcf1, ref, uniqness, threads, outvcf2, bin)
			cmd = "python " + script_path+ "/vcf2json.py" + " -p " + outvcf1 + " -f " + outvcf2 + " -o " + QCVCF
			subprocess.call(cmd, shell=True)

			###report
			run_cmd(
				[
				"python",
				script_path + "/report.py",
				"-f", QCFQ, 
				"-b", QCBAM,
				"-m", QCVCF,
				"-o", htmlfile
				],
				"report",
				)


		line=INFILE.readline()

	INFILE.close()
	
