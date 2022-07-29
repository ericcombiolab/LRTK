import getopt
import multiprocessing
import os
import subprocess
import sys

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

from utility import *

def Phasing_HAPCUT2(bam, vcf, ref, threads, outfile):
	###########################################################################
	# Use HAPCUT2 to phase detected variations                                #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split( outdir + "/" + outfile)
	
	#temp paths and files
	tmp_unlink  = outdir + "/tmp_unlinked_fragments"
	tmp_link    = outdir + "/tmp_linked_fragments"
	
	#Variation calling
	run_cmd(
		[
		"extractHAIRS", 
		"--10x", "1", 
		"--bam", bam, 
		"--VCF", vcf, 
		"--out", tmp_unlink
		],
		"unlink"
	)
	
	run_cmd(
		[
		"python", 
		python_path + "LinkFragments.py",
		"--bam", bam, 
		"--VCF", vcf, 
		"--fragments", tmp_unlink, 
		"--out", tmp_link
		],
		"link"
	)

	run_cmd(
		[
		"HAPCUT2", 
		"--nf", "1", 
		"--fragments", tmp_link, 
		"--VCF", vcf, 
		"--out", outfile
		],
		"hapcut2"
	)	

	#Remove tmp files
	#os.remove(tmp_unlink)
	#os.remove(tmp_link)


def Phasing_WhatsHap(bam, vcf, ref, threads, outfile):
	###########################################################################
	# Use samtools to call the small variations based on the bam files.       #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split( outdir + "/" + outfile)

	###final output
	tmp_vcf     = outdir + "/phased.vcf"

	#Variation calling
	run_cmd(
		[
		"whatshap", 
		"phase", 
		"-o", outfile, 
		"-r", ref, 
		vcf, bam
		],
		"WhatsHap",
	)

	
def Phasing_SpecHap(bam,  vcf, ref, threads, outfile, bin):
	###########################################################################
	# Use GATK to call the small variations based on the bam files.           #
	###########################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split( outdir + "/" + outfile)

	#temp paths and files
	tmp_region  = outdir + "/regions.bed"
	tmp_vcf     = outdir + "/phased.vcf"
	tmp_fai     = ref + ".fai"

	#Variation calling
	#temp paths and files
	tmp_unlink  = outdir + "/tmp_unlinked_fragments"
	tmp_link    = outdir + "/tmp_linked_fragments"
	sort_link   = outdir + "/tmp_sort_fragments"
	bcbed       = outdir + "/barcode_spnanning.bed"
	bcbedz      = outdir + "/barcode_spnanning.bed.gz"
	
	#
	cmd = "chmod 755 " + python_path + "SpecHap/SpecHap"
	subprocess.call(cmd, shell=True)
	cmd = "chmod 755 " + python_path + "SpecHap/BarcodeExtract"
	subprocess.call(cmd, shell=True)

	#Variation calling
	subprocess.call(["extractHAIRS", "--10x", "1", "--bam", bam, "--VCF", vcf, "--out", tmp_unlink])
	subprocess.call(["sort -n -k6 ", tmp_unlink  " > ", sort_link])
	subprocess.call([python_path + "SpecHap/BarcodeExtract", bam, bcbed])
	subprocess.call(["bgzip -c ", bcbed, " > ", bcbedz])
	subprocess.call(["tabix -p bed ",bcbedz])
	subprocess.call([python_path + "SpecHap/SpecHap", "--vcf", vcf, "--frag", tmp_unlink, "--out", tmp_vcf, "--10x"])

	#Remove tmp files
	#os.remove(tmp_unlink)
	#os.remove(tmp_link)

