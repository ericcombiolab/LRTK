import getopt
import multiprocessing
import os
import subprocess
import sys
from utility import *

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
	cmd="extractHAIRS" + " --10x " + str("1") + " --bam " + bam + " --VCF " + vcf + " --out " + tmp_unlink
	run_cmd(cmd,"extractHAIRS")

	cmd="python " + python_path + "LinkFragments.py" + " --bam " + bam + " --VCF " + vcf + " --fragments " + tmp_unlink + " --out " + tmp_link
	run_cmd(cmd,"LinkFragments")

	cmd="HAPCUT2 " + " --nf " + str("1") + " --fragments " + tmp_link + " --VCF " + vcf + " --out " + outfile
	run_cmd(cmd,"hapcut2")

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
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	###final output
	tmp_vcf = outdir + "/phased.vcf"

	#Variation calling
	cmd="whatshap" + " phase " + " -o " + outfile + " -r " + ref + " " + vcf + " " + bam
	run_cmd(cmd,"WhatsHap")

def Phasing_WhatsHap_meta(bam, vcf, ref, ploidy, threads, outfile):
	###########################################################################
	# Use whatshap to phase heterozygous SNVs                                 #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split( outdir + "/" + outfile)

	cmd="whatshap polyphase " + vcf + " " + bam + " --ploidy " + ploidy + " --reference " + ref + " -o " + outfile
	run_cmd(cmd, "polyphase")

