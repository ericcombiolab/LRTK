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
import metagenome_assembly
from utility import *


def moduleFQ(args):	
	# To convert different fastq formats
	# convert 10x to ULRF
	if(args.input_type == "10x"): #convert 10x to ULRF
		conversion.TENx2ULRF(args.input_fastq1, args.input_fastq2, args.output_fastq1, args.output_fastq2, args.barcodes, args.host, args.filter, args.sort, args.threads, args.genome)
	elif(args.input_type == "stLFR"): # convert stLFR to ULRF
		conversion.stLFR2ULRF(args.input_fastq1, args.input_fastq2, args.output_fastq1, args.output_fastq2, args.barcodes, args.host, args.filter, args.sort, args.threads, args.genome)
	elif(args.input_type == "TELLSeq"): # convert TELL-Seq to ULRF
		conversion.TELLSeq2ULRF(args.input_fastq1, args.input_fastq2, args.index_fastq, args.output_fastq1, args.output_fastq2, args.host, args.filter, args.sort, args.threads, args.genome)
	else:
		print("Error: unsupported format conversion\n")

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
	if((args.platform == "10x") and (args.genome == "human")):
		alignment.align_10x(args.input_fastq1, args.input_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.threads)
	elif(args.platform == "stLFR" and (args.genome == "human")):
		alignment.align_stLFR(args.input_fastq1, args.input_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.threads)
	elif(args.platform == "TELLSeq" and (args.genome == "human")):
		alignment.align_TELLSeq(args.input_fastq1, args.input_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.threads)
	elif(args.genome == "metagenome"):
		alignment.metagenome_align(args.input_fastq1, args.input_fastq2, args.read_group, args.reference, args.outfile, args.sort, args.mark_duplication, args.platform, args.threads)
	else:
		print("Error: unknown linked-read technology!")


def moduleSNV(args):
	if(args.genome == "human"):
		if((args.application == "FreeBayes")):
			small_variants.calling_freebayes(args.bam, args.reference, args.threads, args.outfile)
		elif(args.application == "Samtools"):
			small_variants.calling_samtools(args.bam, args.reference, args.threads, args.outfile)
		elif(args.application == "GATK"):
			small_variants.calling_GATK(args.bam, args.reference, args.threads, args.outfile)
		else:
			print("Error: unknown SNV/INDEL caller!")
	elif(args.genome == "metagenome"):
		if((args.application == "FreeBayes")):
			small_variants.calling_metasnv_freebayes(args.bam, args.reference, args.threads, args.outfile)
		elif(args.application == "Samtools"):
			small_variants.calling_metasnv_samtools(args.bam, args.reference, args.threads, args.outfile)
		elif(args.application == "instrain"):
			small_variants.calling_metasnv_inStrain(args.bam, args.reference, args.threads, args.outfile)
		else:
			print("Error: unknown SNV/INDEL caller!")
	else:
		print("Error: unknown organism!")


def moduleSV(args):
	if(args.application == "Aquila"):
		print("Here is moduleSV")
		large_variants.calling_Aquila(args.bam, args.vcf, args.reference, args.uniqness, args.threads, args.outfile)
	elif(args.application == "LinkedSV"):
		large_variants.calling_LinkedSV(args.bam, args.reference, args.threads, args.outfile)
	elif(args.application == "VALOR"):
		large_variants.calling_VALOR(args.bam, args.reference, args.sonic, args.threads, args.outfile)
	else:
		print("Error: unknown SV caller!")

def modulePHASE(args):
	###make preparations
	if(args.genome == "human"):
		if((args.application == "HapCUT2")):
			phasing.Phasing_HAPCUT2(args.bam, args.vcf, args.reference, args.threads, args.outfile)
		elif(args.application == "WhatsHap"):
			phasing.Phasing_WhatsHap(args.bam, args.vcf, args.reference, args.threads, args.outfile)
		else:
			print("Error: unknown phasing tool for human genome!")    
	elif(args.genome == "metagenome"):
		if((args.application == "WhatsHap")):
			phasing.Phasing_WhatsHap_meta(args.bam, args.vcf, args.reference, args.number, args.threads, args.outfile)
		else:
			print("Error: unknown phasing tool for metagenome!")   
	else:
		print("Error: unknown organism!")

def moduleASSEMBLY(args):
	metagenome_assembly.MAGconstruction(args.fq1, args.fq2, args.metaspades, args.athena_l, args.athena_h, args.lt, args.outfile, args.threads)

def moduleWGS(args):
	###FASTQ processing && mapping
	sample_infile = args.sample_info
	rg = args.read_group
	INFILE = open(sample_infile,"rt")
	line = INFILE.readline()
	humanref=args.database + "/GRCH38/genome.fa"
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
				conversion.TENx2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, humanref, "Yes", "Yes", threads, "human")
				alignment.align_10x(outfq1, outfq2, rg, ref, bamfile, "Yes", "Yes", threads)
				
			elif(tech == "stLFR"):
				whitelist=args.database + "/WhiteList/white_list_stlfr_barcode.fa"
				conversion.stLFR2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, humanref, "Yes", "Yes", threads, "human")
				alignment.align_stLFR(outfq1, outfq2, rg, ref, bamfile, "Yes", "Yes", threads)

			elif(tech == "TELLSeq"):
				conversion.TELLSeq2ULRF(input_fastq1, input_fastq2, index_fastq, outfq1, outfq2, humanref,"Yes", "Yes", threads,"human")
				alignment.align_TELLSeq(outfq1, outfq2, rg, ref, bamfile, "Yes", "Yes", threads)

			else:
				print("Unknown linked-read technology!")

			samtools_path = shutil.which("samtools")
			cmd = "chmod 755 " + script_path + "long_fragment/construct_fragment"
			subprocess.call(cmd, shell=True)
			
			cmd = script_path + "long_fragment/construct_fragment" + " -d " + str(200000) + " -t " +  str(threads) + " --bam " + bamfile + " --output " + resBAM + "/constructedfragment.xls" + " 1>" + resBAM + "/constructedfragment.stat.xls"
			run_cmd(cmd, "construct_fragment")
			cmd = "python " + script_path+ "/Fragment2json.v2.py" + " -b " + bamfile + " -f " + resBAM + "/constructedfragment.xls" + " -i " + bamfile + ".insert.list" + " -s " + resBAM + "/constructedfragment.stat.xls" + " -t " + samtools_path + " -o " + QCBAM
			subprocess.call(cmd, shell=True)

			###Variants calling && phasing
			small_variants.calling_freebayes(bamfile, ref, threads, outvcf1)
			phasing.Phasing_HAPCUT2(bamfile, outvcf1, ref, threads, outvcf3)
			large_variants.calling_Aquila(bamfile, outvcf1, ref, uniqness, threads, outvcf2)
			cmd = "python " + script_path+ "/vcf2json.py" + " -p " + outvcf1 + " -f " + outvcf2 + " -o " + QCVCF
			subprocess.call(cmd, shell=True)

			###report
			cmd="python" + " " + script_path + "/report.py" + " -f " + QCFQ + " -b " + QCBAM + " -m " + QCVCF + " -o " + htmlfile
			subprocess.call(cmd, shell=True)

		line=INFILE.readline()

	INFILE.close()

def moduleMWGS(args):
	###FASTQ processing && mapping

	bin = script_path + "/.."

	sample_infile = args.sample_info
	multi_infile = args.multi_info
	rg=args.read_group
	threads=args.threads
	humanref=args.database + "/GRCH38/genome.fa"
	metaref=args.database + "/UHGG/UHGG.reference.fa"

	outdir=args.outdir
	outdir_multi = outdir + "/" + "Merge"
	outdir_multi_a = outdir + "/Merge/Analysis"
	outdir_multi_r = outdir + "/Merge/Report"

	INFILE = open(sample_infile,"rt")
	line = INFILE.readline()

	###processing single sample
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
			resPS=outdir_sample+"/"+"Phasing"
			resRP=outdir_sample+"/"+"Report"

			os.makedirs(outdir_sample)
			os.makedirs(resFQ)
			os.makedirs(resBAM)
			os.makedirs(resQC)
			os.makedirs(resSNV)
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
			QCFQ   = resQC + "/" + "FASTQ.QC.json"

			bamfile = resBAM + "/" + "align.round2.bam"
			bamstat = outdirbam + ".Clean.sort.rmdup.bam.AlignmentStat.xls"
			QCBAM   = resQC + "/" + "Fragment.json"

			outvcf1 = resSNV + "/" + sample_name+".smallvariants.vcf"
			outvcf3 = resPS  + "/" + sample_name+".phased.vcf"
			QCVCF   = resQC  + "/" + "VCF.json"

			iter_ref = resBAM + "/" + "UHGG.reference.subset.fa"
			htmlfile = resRP + "/" + sample_name + ".report.html"

			if((tech == "10x")):
				whitelist=args.database + "/WhiteList/white_list_10x_barcode.fa"
				conversion.TENx2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, humanref, "Yes", "Yes", threads, "metagenome")
				alignment.metagenome_align(outfq1, outfq2, rg, metaref, bamfile, "Yes", "Yes","10x", threads)

			elif(tech == "stLFR"):
				whitelist=args.database + "/WhiteList/white_list_stlfr_barcode.fa"
				conversion.stLFR2ULRF(input_fastq1, input_fastq2, outfq1, outfq2, whitelist, humanref, "Yes", "Yes", threads, "metagenome")
				alignment.metagenome_align(outfq1, outfq2, rg, metaref, bamfile, "Yes", "Yes", "stlfr",threads)

			elif(tech == "TELLSeq"):
				conversion.TELLSeq2ULRF(input_fastq1, input_fastq2, index_fastq, outfq1, outfq2, humanref, "Yes", "Yes", threads, "metagenome")
				alignment.metagenome_align(outfq1, outfq2, rg, metaref, bamfile, "Yes", "Yes", "tellseq", threads)

			else:
				print("Unknown linked-read technology!")

			###QC module
			cmd = "cp " + resFQ + "/FASTQ.QC.json " + QCFQ
			subprocess.call(cmd, shell=True)
			samtools_path = shutil.which("samtools")
			
			cmd = "chmod 755 " + script_path + "long_fragment/construct_fragment"
			subprocess.call(cmd, shell=True)
			
			cmd = script_path + "long_fragment/construct_fragment" + " -d " + str(200000) + " -t " +  str(threads) + " --bam " + bamfile + " --output " + resBAM + "/constructedfragment.xls" + " 1>" + resBAM + "/constructedfragment.stat.xls"
			run_cmd(cmd, "construct_fragment")
			
			cmd = "python " + script_path+ "/Fragment2json.v2.py" + " -b " + bamfile + " -f " + resBAM + "/constructedfragment.xls" + " -i " + resBAM + "/constructedfragment.xls.insert_size" + " -s " + resBAM + "/constructedfragment.stat.xls" + " -t " + samtools_path + " -o " + QCBAM
			subprocess.call(cmd, shell=True)

			###Variants calling && phasing
			small_variants.calling_metasnv_samtools(bamfile, iter_ref, threads, outvcf1)
			cmd = "python " + script_path+ "/MetaSNV2json.py" + " -a " + resBAM + "/UHGG.coverage_average.iter.xls" + " -f " + resBAM + "/UHGG.reference.subset.fa.fai" + " -v " + outvcf1 + " -o " + QCVCF
			subprocess.call(cmd, shell=True)

			###report
			cmd = "python " + script_path + "/MetaReport.py" + " -f " + QCFQ + " -b " + QCBAM + " -m " + QCVCF + " -o " + htmlfile
			subprocess.call(cmd, shell=True)

		line = INFILE.readline()

	INFILE.close()
	###processing multiple samples
	if multi_infile is not None:
		#print(multi_infile)
		
		cmd = "python " + script_path+ "/MultipleMetagenomicAnalysis.py -i " + multi_infile + " -a " + outdir + " -d " + args.database + " -o " + outdir_multi
		#print(cmd)
		subprocess.call(cmd, shell=True)

		cmd = "python " + script_path+ "/MultipleReport.py -a " + outdir_multi + "/report/ABD/Abundance.stastistics.xls" + " -s " + outdir_multi + "/report/SNV/SNV.statistics.xls" +  " -o " + outdir_multi + "/report/multi.html"
		#print(cmd)
		subprocess.call(cmd, shell=True)

