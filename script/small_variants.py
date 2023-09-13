import getopt
import multiprocessing
import os
import subprocess
import sys

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

from utility import *

def calling_freebayes(bam, ref, threads, outfile):
	###########################################################################
	# Use freebayes to call the small variations.                             #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	#temp paths and files
	tmp_region  = outdir + "/regions.bed"
	tmp_vcf     = outdir + "/tmp.vcf"

	#Variation calling
	#logging("calling SNP/INDEL start.")
	cmd="python " + python_path + "fasta_generate_regions.py " + ref + " 100000 " + tmp_region
	run_cmd(cmd,"generate_bed")

	cmd="freebayes-parallel " + tmp_region + " " + str(threads) + " -f " + ref + " " + bam + " > " + tmp_vcf
	run_cmd(cmd, "SNV_freebayes")

	#Filteration
	cmd="perl " +python_path + "Filt.SmallVariants.pl " + tmp_vcf + " " + outfile
	subprocess.call(cmd, shell=True)

	os.remove(tmp_vcf)
	logging("calling SNP/INDEL end.")

def calling_samtools(bam, ref, threads, outfile):
	###########################################################################
	# Use samtools to call the small variations.                              #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	#temp paths and files
	tmp_bcf     = outdir + "/tmp.bcf"
	tmp_vcf     = outdir + "/tmp.vcf"

	#Variation calling
	cmd = "bcftools " + "mpileup " + "-f " + ref + " -q 20 -C 50 -Q 13 -E --threads " + str(threads) + " " + bam  + " | bcftools call -cv -Ob --threads " + str(threads)  + " -o " + tmp_bcf
	run_cmd(cmd, "SNV_samtools")

	#Filteration
	cmd = "bcftools filter -sLowQual -g3 -G10 -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15)' " + tmp_bcf + " -o " + tmp_vcf
	run_cmd(cmd, "SNV_filtering")
	subprocess.call(["perl", python_path + "Filt.SmallVariants.pl", "", "-i", tmp_vcf, "-o", outfile])

	
def calling_GATK(bam, ref, threads, outfile):
	###########################################################################
	# Use GATK to call the small variations based on the bam files.           #
	###########################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	#temp paths and files
	tmp_vcf     = outdir + "/tmp.vcf"

	#Variation calling
	cmd="gatk3 " + " -T " + " UnifiedGenotyper " + " -R " + ref + " -I " + bam + "-o" + tmp_vcf + "-stand_call_conf 50 -A RMSMappingQuality -baq CALCULATE_AS_NECESSARY"
	run_cmd(cmd, "SNV_GATK")

	#Filteration
	#cmd="gatk3" + " -T " + "VariantFiltration" + " -R " + ref + " -V " + tmp_vcf + " -o " + outfile + " --filterExpression " + "'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)'" + " --filterName " + "HARD_TO_VALIDATE" + " --filterExpression " + "'DP < 3'" + " --filterName " + "LOW_READ_SUPPORT"
	cmd="gatk3" + " -T " + "VariantFiltration" + " -R " + ref + " -V " + tmp_vcf + " -o " + outfile  + " --filterExpression " + "'DP < 3'" + " --filterName " + "LOW_READ_SUPPORT"
	subprocess.call(cmd, "SNV_filtering")

def calling_metasnv_inStrain(bam, ref, threads, outfile):
	###########################################################################
	#Use inStrain to call the SNVs based on the bam files.                    #
	###########################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	cmd="inStrain profile " + bam + " " + ref + " -o " + outfile +  " -p " + str(threads)
	run_cmd(cmd, "metaSNV_inStrain")

def calling_metasnv_freebayes(bam, ref, threads, outfile):
	###########################################################################
	#Use freebayes to call the SNVs based on the bam files.                    #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	cmd="freebayes -f " + ref + " " + bam + " > " + outfile
	run_cmd(cmd, "metaSNV_freebayes")

def calling_metasnv_samtools(bam, ref, threads, outfile):
	###########################################################################
	#Use samtools to call the SNVs based on the bam files.                   #
	###########################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	tmp_bcf = outdir + "/tmp.bcf"

	cmd="bcftools mpileup -f " + ref + " -q 20 -C 50 -Q 13 -d 1000000 -I --threads  " + str(threads) + " " + bam + " |bcftools call -cv -Ou --threads " + str(threads) + " -o " + tmp_bcf
	run_cmd(cmd, "metaSNV_samtools")
	cmd="bcftools  filter -sLowQual -g3 -G10 " + tmp_bcf + " -o " + outfile
	run_cmd(cmd, "SNV_filtering")

	os.remove(tmp_bcf)

