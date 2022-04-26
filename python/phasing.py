import getopt
import multiprocessing
import os
import subprocess

def Phasing_HAPCUT2(bam, vcf, ref, threads, outfile, bin):
	###########################################################################
	# Use HAPCUT2 to phase detected variations                                #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
	
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
		bin + "/python/LinkFragments.py",
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


def Phasing_WhatsHap(bam, vcf, ref, threads, outfile, bin):
	###########################################################################
	# Use samtools to call the small variations based on the bam files.       #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()

	###final output
	tmp_vcf     = outdir + "/phased.vcf"

	#Variation calling
	run_cmd(
		[
		"whatshap", 
		"phase", 
		"-o", tmp_vcf, 
		"--reference=", ref, 
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

	#temp paths and files
	tmp_region  = outdir + "/regions.bed"
	tmp_vcf     = outdir + "/tmp.vcf"
	tmp_fai     = ref + ".fai"

	#Variation calling
	#temp paths and files
	tmp_unlink  = outdir + "/tmp_unlinked_fragments"
	tmp_link    = outdir + "/tmp_linked_fragments"

	#Variation calling
	subprocess.call(["extractHAIRS", "--10x", "1", "--bam", bam, "--VCF", vcf, "--out", tmp_unlink])
	subprocess.call(["python", bin + "/SCRIPTS/LinkFragments.py", "--bam", bam, "--VCF", vcf, "--fragments", tmp_unlink, "--out", tmp_link])
	subprocess.call(["HAPCUT2", "--nf", "1", "--fragments", tmp_link, "--VCF", vcf, "--out", outfile])

	#Remove tmp files
	#os.remove(tmp_unlink)
	#os.remove(tmp_link)

