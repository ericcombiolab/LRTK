import getopt
import multiprocessing
import os
import subprocess
from utility import *

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

def calling_Aquila(bam, vcf, ref, umap, threads, outfile):
	###########################################################################
	# Use Aquila to call the large structural variations.                     #
	###########################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	#temp paths and files
	tmp_region  = outdir + "/regions.bed"
	tmp_vcf     = outdir + "/tmp.vcf"
	tmp_fai     = ref + ".fai"

	for i in range(23):
		chrom= str(i+1)
		outdirchr=outdir + "/CHR" + chrom

		cmd="Aquila_step1" + " --bam_file " + bam + " --vcf_file " + vcf + " --sample_name " + "Aquila " + "--chr_start " + chrom + " --chr_end " + chrom + " --uniq_map_dir " + umap + " --num_threads_for_samtools_sort " + str(threads) + " --out_dir " + outdirchr + " --num_threads " + str("2")
		run_cmd(cmd,"SV_Aquila_step1")

		cmd="Aquila_step2" + " --reference " + ref + " --num_threads " + str(threads) + " --chr_start " + chrom + " --chr_end " + chrom + " --out_dir " + outdirchr + " --num_threads_spades " + str("2")
		run_cmd(cmd,"SV_Aquila_step2")

		cmd="Aquila_assembly_based_variants_call" + " --assembly_dir " + outdirchr + " --out_dir " + outdirchr + "/VariantsResults" + " --ref_file " + ref + " --num_of_threads " + str(threads) + " --chr_start " + chrom + " --chr_end " + chrom + " --all_regions_flag " + str("1")
		run_cmd(cmd,"SV_Aquila_step3")

		cmd="Aquila_clean" + " --assembly_dir " + outdirchr
		run_cmd(cmd,"SV_Aquila_clean")

	##Merge final vcf
	cmd="cat " 
	for i in range(23):
		chrom= str(i+1)
		cmd= cmd + outdir + "/CHR" + chrom + "/VariantsResults/Aquila_final_sorted.vcf "

	cmd=cmd+ " > " + tmp_vcf
	subprocess.call(cmd, shell=True)
	cmd="python " + python_path + "Reformat.py" + " -r " + ref + " -i " + tmp_vcf + " -o " + outfile + " --add_header 38 --base_norm --gz_tbi"
	subprocess.call(cmd, shell=True)


def calling_LinkedSV(bam, ref, threads, outfile):
	######################################################################################
	# Use LinkedSV to call the large structural variations based on the bam files.       #
	######################################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	linkedsv= python_path + "LinkedSV/linkedsv.py"
	#Variation calling
	cmd="python " + linkedsv + " -i " + bam + " -d " + outdir + " -r " + ref + " -v " + " hg38 " + " -t " + str(threads) + " --germline_mode"
	run_cmd(cmd, "SV_LinkedSV")

	
def calling_VALOR(bam, ref, sonic, threads, outfile):
	#######################################################################################
	# Use VALOR to call the large structural variations based on the bam files.           #
	#######################################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	valor = python_path + "valor/valor"
	cmd = "chmod 755 " + valor
	subprocess.call(cmd, shell=True)
	#Variation calling
	subprocess.call([valor, "-i", bam, "-o", outdir, "-s", sonic, "-f", "INV,DUP,IDUP,TRA,ITRA,DEL", "-t", threads])

def calling_LGT(genomes, anno, threads, outfile):
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)

	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	cmd="MetaCHIP PI -p " +  outdir + " -r c -t " + str(threads) + " -i " + SEQDIR + " -x fasta -taxon " + SEQANNO
	subprocess.call(cmd, shell=True)
	cmd="MetaCHIP BP -p " + outdir + " -r c -t " + str(threads)
	subprocess.call(cmd, shell=True)

