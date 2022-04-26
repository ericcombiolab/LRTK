import getopt
import multiprocessing
import os
import subprocess

def calling_Aquila(bam, vcf, ref, umap, threads, outfile, bin):
	###########################################################################
	# Use Aquila to call the large structural variations.                     #
	###########################################################################

	print("Here is Large Variants module")
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()

	#temp paths and files
	tmp_region  = outdir + "/regions.bed"
	tmp_vcf     = outdir + "/tmp.vcf"
	tmp_fai     = ref + ".fai"

	for i in range(23):
		chrom= str(i+1)
		outdirchr=outdir + "/CHR" + chrom
		
		run_cmd(
			[
			"Aquila_step1", 
			"--bam_file", bam, 
			"--vcf_file", vcf, 
			"--sample_name", "Aquila", 
			"--chr_start", chrom, 
			"--chr_end", chrom, 
			"--uniq_map_dir", umap, 
			"--num_threads_for_samtools_sort", threads, 
			"--out_dir", outdirchr, 
			"--num_threads", "5"
			],
			"step1"
		)
		
		run_cmd(
			[
			"Aquila_step2", 
			"--reference", ref, 
			"--num_threads", "5", 
			"--chr_start", chrom, 
			"--chr_end", chrom, 
			"--out_dir", outdirchr, 
			"--num_threads_spades", "5"
			],
			"step2"
		)
		
		run_cmd(
			[
			"Aquila_assembly_based_variants_call", 
			"--assembly_dir", outdirchr, 
			"--out_dir", outdirchr + "/VariantsResults", 
			"--ref_file", ref, 
			"--num_of_threads", threads, 
			"--chr_start", chrom, 
			"--chr_end", chrom, 
			"--all_regions_flag", "1"
			],
			"variants_call"
		)

		run_cmd(
			[
			"Aquila_phasing_all_variants", 
			"--assembly_vcf", outdirchr + "/VariantsResults/Aquila_final_sorted.vcf", 
			"--vcf_file", vcf, 
			"--assembly_dir", outdirchr,
			"--chr_start", chrom, 
			"--chr_end", chrom, 
			"--out_dir", outdirchr + "/Aquila_Phasing_Results"
			],
			"phasing_variants"
		)
	
		run_cmd(
			[
			"Aquila_clean", 
			"--assembly_dir", outdirchr
			],
			"clean"
		)


def calling_LinkedSV(bam, ref, threads, outfile, bin):
	######################################################################################
	# Use LinkedSV to call the large structural variations based on the bam files.       #
	######################################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()

	#Variation calling
	subprocess.call(["python", "linkedsv.py", "-i", bam, "-d", outdir, "-r", ref, "-v", "hg38", "-t", "64", "--germline_mode"])
	#subprocess.call(["freebayes-parallel", tmp_region, threads, "-f", ref, bam, ">", tmp_vcf])

	#Filteration
	#subprocess.call(["samtools", "aln", white_list, bc_fq, "-f", bc_sai])

	
def calling_VALOR(bam, ref, sonic, threads, outfile, bin):
	#######################################################################################
	# Use VALOR to call the large structural variations based on the bam files.           #
	#######################################################################################

	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()


	#Variation calling
	subprocess.call(["valor", "-i", bam, "-o", outdir, "-s", sonic, "-f", "DUP,DEL", "-t", threads])
	#subprocess.call(["freebayes-parallel", tmp_region, threads, "-f", ref, bam, ">", tmp_vcf])

	#Filteration
	#subprocess.call(["samtools", "aln", white_list, bc_fq, "-f", bc_sai])	
