import os
import sys
import subprocess

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

from utility import *

def align_10x(bq1, bq2, fq1, fq2, rg, ref, outfile, sort, mark, threads, bin):
	###########################################################################
	# To align the clean ULRF sequencing reads to reference genome            #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
	
	#temp paths and files
	tmp_ema_sam  = outdir + "/tmp.ema.sam"
	tmp_ema_bam  = outdir + "/tmp.ema.bam"
	tmp_bwa_sam  = outdir + "/tmp.bwa.sam"
	tmp_bwa_bam  = outdir + "/tmp.bwa.bam"
	tmp_merge_bam = outdir + "/tmp.merge.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"

	rg="@RG\\tID:example\\tSM:example"

	#Align reads
	pipe = subprocess.Popen(
			[
			python_path+"../src/EMA/ema", 
			"align",
			"-1", bq1,
			"-2", bq2,
			"-r", ref,
			"-R", rg,
			"-p", "10x",
			"-t", threads
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.DEVNULL,
		)

	subprocess.check_output(
			[
			"samtools",
			"view",
			"-Sb",
			"-",
			"-o",
			tmp_ema_bam,
			],
			stdin=pipe.stdout,
		)
	pipe.communicate()

	pipe = subprocess.Popen(
			[
			"bwa", 
			"mem", 
			ref, 
			fq1, 
			fq2, 
			"-R", rg, 
			"-t", threads,
			], 
			stdout=subprocess.PIPE,	
			stderr=subprocess.DEVNULL,
		)

	subprocess.call(
			[
			"samtools", 
			"view", 
			"-Sb", 
			"-",
			"-o", tmp_bwa_bam
			],
			stdin=pipe.stdout,
		)
	pipe.communicate()

	#Bam Post Processing
	run_cmd(
		[
		"samtools",
		"merge",
		tmp_merge_bam, 
		tmp_bwa_bam, 
		tmp_ema_bam
		],
		"merge_bam",
	)
	
	run_cmd(
		[
		"samtools", 
		"sort", 
		"--threads", threads, 
		tmp_merge_bam, 
		"-o", tmp_sort_bam
		],
		"sort_bam",
	)

	run_cmd(
		[
		"picard", 
		"MarkDuplicates", 
		"I=", tmp_sort_bam, 
		"O=", outfile, 
		"M=", tmp_markdup_mat, 
		"BARCODE_TAG=","BC"
		],
		"mark_dup_bam",
	)

	run_cmd(
		[
		"samtools", 
		"index", 
		outfile
		],
		"index_bam"
	)	

	#QC
	run_cmd(
		[
		python_path + "../src/long_fragment/construct_fragment",
		"--bam", outfile,
		"--output", outdir
		],
		"fragment_stat",
	)
	
	run_cmd(
		[
		"perl", python_path+"AlignStat_WGS.pl", 
		"-in", outfile, 
		"-sam", filename, 
		"-samtools", "/home/comp/cschaoyang/SOFTWARE/ANACONDA3/V2019/bin/samtools", 
		"-outdir", outdir, 
		"-move", outdir
		],
		"bam_stat"
	)
	
	#remove tmp files
	os.remove(tmp_ema_bam)
	os.remove(tmp_bwa_bam)
	os.remove(tmp_merge_bam)
	os.remove(tmp_sort_bam)


def align_stLFR(bq1, bq2, fq1, fq2, rg, ref, outfile, sort, mark, threads, bin):
	###########################################################################
	# To align clean ULRF (bc=30bp) sequencing reads to reference genome      #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()

	#temp paths and files
	tmp_ema_sam  = outdir + "/tmp.ema.sam"
	tmp_ema_bam  = outdir + "/tmp.ema.bam"
	tmp_bwa_sam  = outdir + "/tmp.bwa.sam"
	tmp_bwa_bam  = outdir + "/tmp.bwa.bam"
	tmp_merge_bam = outdir + "/tmp.merge.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"

	rg="@RG\\tID:example\\tSM:example"
	
	#Align reads
	pipe = subprocess.Popen(
		[
		python_path+"../src/EMA/ema",
		"align",
		"-1", bq1,
		"-2", bq2,
		"-r", ref, 
		"-R", rg,
		"-p", "stLFR",
		"-t", threads,
		],
		stdout=subprocess.PIPE,
		stderr=subprocess.DEVNULL,
	)
	
	subprocess.check_output(
		[
		"samtools",
		"view",
		"-Sb",
		"-",
		"-o",
		tmp_ema_bam,
		],
		stdin=pipe.stdout,
	)
	
	pipe.communicate()

	pipe = subprocess.Popen(
			[
			"bwa",
			"mem",
			ref,
			fq1,
			fq2,
			"-R", rg,
			"-t", threads,
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.DEVNULL,
		)

	subprocess.call(
			[
			"samtools",
			"view",
			"-Sb",
			"-",
			"-o", tmp_bwa_bam,
			],
			stdin=pipe.stdout,
		)

	pipe.communicate()

	#Bam Post Processing
	run_cmd(
		[
		"samtools",
		"merge",
		tmp_merge_bam,
		tmp_bwa_bam,	
		tmp_ema_bam
		],
		"merge_bam",
	)

	run_cmd(
		[
		"samtools",
		"sort",
		"--threads", threads,
		tmp_merge_bam,
		"-o", tmp_sort_bam,
		],
		"sort_bam",
	)

	run_cmd(
		[
		"picard",
		"MarkDuplicates",
		"I=",tmp_sort_bam,
		"O=",outfile,
		"M=",tmp_markdup_mat,
		"BARCODE_TAG=","BC"
		],
		"mark_dup_bam",
	)

	run_cmd(
		[
		"samtools",
		"index",
		outfile
		],
		"index_bam"
	)
	
	#QC
	run_cmd(
		[
		python_path + "../src/long_fragment/construct_fragment",
		"--bam", outfile,
		"--output", outdir,
		],
		"fragment_stat",
	)
	
	run_cmd(
		[
		"perl", python_path+"AlignStat_WGS.pl",
		"-in", outfile,
		"-sam", filename,
		"-samtools","/home/comp/cschaoyang/SOFTWARE/ANACONDA3/V2019/bin/samtools",
		"-outdir", outdir, 
		"-move", outdir,
		],
		"bam_stat"
	)

	#remove tmp files
	os.remove(tmp_ema_sam)
	os.remove(tmp_bwa_bam)
	os.remove(tmp_ema_bam)
	os.remove(tmp_merge_bam)
	os.remove(tmp_sort_bam)


def align_TELLSeq(bq1, bq2, rg, ref, outfile, sort, mark, threads, bin):
	###########################################################################
	# To align clean ULRF (bc=18bp) sequencing reads to reference genome      #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()

	#temp paths and files
	tmp_ema_sam  = outdir + "/tmp.ema.sam"
	tmp_ema_bam  = outdir + "/tmp.ema.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"
	rg="@RG\\tID:example\\tSM:example"

	#Align reads
	pipe = subprocess.Popen(
			[
			python_path+"../src/EMA/ema",
			"align",
			"-1", bq1,
			"-2", bq2,
			"-r", ref,
			"-R", rg,
			"-p", "TELLSeq",
			"-t", threads
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.DEVNULL,
		)

	subprocess.check_output(
		[
		"samtools",
		"view",
		"-Sb",
		"-",
		"-o",
		tmp_ema_bam,
		],
		stdin=pipe.stdout,
	)

	pipe.communicate()
	
	#Bam Post Processing
	run_cmd(
		[
		"samtools",
		"merge",
		tmp_merge_bam,
		tmp_bwa_bam,
		tmp_ema_bam
		],
		"merge_bam",
	)

	run_cmd(
		[
		"samtools",
		"sort",
		"--threads", threads,
		tmp_merge_bam,
		"-o", tmp_sort_bam
		],
		"sort_bam"
	)

	run_cmd(
		[
		"picard",
		"MarkDuplicates",
		"I=", tmp_sort_bam,
		"O=", outfile,
		"M=", tmp_markdup_mat,
		"BARCODE_TAG=","BC",
		],
		"mark_dup_bam",
	)

	run_cmd(
		[
		"samtools",
		"index",
		outfile
		],
		"index_bam"
	)

	#QC
	run_cmd(
		[
		python_path + "../src/long_fragment/construct_fragment",
		"--bam", outfile,
		"--output", outdir,
		],
		"fragment_stat",
	)
	
	run_cmd(
		[
		"perl", python_path+"AlignStat_WGS.pl",
		"-in", outfile,
		"-sam", filename,
		"-samtools","/home/comp/cschaoyang/SOFTWARE/ANACONDA3/V2019/bin/samtools",  
		"-outdir", outdir,
		"-move", outdir
		],
		"bam_stat",
	)

	#remove tmp files
	os.remove(tmp_ema_sam)
	os.remove(tmp_ema_bam)
	os.remove(tmp_sort_bam)

