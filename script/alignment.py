import os
import sys
import subprocess
import shutil

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

from utility import *

def align_10x(fq1, fq2, rg, ref, outfile, sort, mark, threads):
	###########################################################################
	# To align the clean ULRF sequencing reads to reference genome            #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir,filename)=os.path.split(outdir + "/" + outfile)
	
	#temp paths and files
	tmp_b_fq1 = outdir + "/tmp.barcoded.1.fq"
	tmp_b_fq2 = outdir + "/tmp.barcoded.2.fq"
	tmp_u_fq1 = outdir + "/tmp.unbarcoded.1.fq"
	tmp_u_fq2 = outdir + "/tmp.unbarcoded.2.fq"
	tmp_ema_sam  = outdir + "/tmp.ema.sam"
	tmp_ema_bam  = outdir + "/tmp.ema.bam"
	tmp_bwa_sam  = outdir + "/tmp.bwa.sam"
	tmp_bwa_bam  = outdir + "/tmp.bwa.bam"
	tmp_merge_bam = outdir + "/tmp.merge.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"
	rg="\'"+rg+"\'"

	#split and sort reads
	cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq1 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_b_fq1
	cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq2 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_b_fq2
	cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq1 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_u_fq1
	cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq2 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_u_fq2
	run_cmd(cmd1, "sort_barcoded_FQ1")
	run_cmd(cmd2, "sort_unbarcoded_FQ1")
	run_cmd(cmd3, "sort_barcoded_FQ2")
	run_cmd(cmd4, "sort_unbarcoded_FQ2")

	#Align reads
	cmd = "chmod 755 " + python_path + "EMA/ema"
	subprocess.call(cmd, shell=True)
	cmd = python_path + "EMA/ema" + " align "+" -1 "+ tmp_b_fq1 + " -2 "+ tmp_b_fq2 +" -r "+ ref + " -R "+ rg + " -p "+ " 10x "+" -t "+ str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_ema_bam + " 2>/dev/null"
	run_cmd(cmd, "ema_alignment")
	cmd = "bwa " + " mem " + ref + " " + tmp_u_fq1 + " " + tmp_u_fq2 + " -R " + rg + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_bwa_bam + " 2>/dev/null"
	run_cmd(cmd, "bwa_alignment")

	#Bam Post Processing
	cmd = "samtools " + " merge " + " -f " + tmp_merge_bam + " " + tmp_bwa_bam + " " + tmp_ema_bam + " --threads " + str(threads)
	run_cmd(cmd, "merge_bam")

	cmd="samtools " + " sort " + tmp_merge_bam + " -o " + tmp_sort_bam + " --threads " + str(threads)
	run_cmd(cmd, "sort_bam")

	cmd="picard " + " MarkDuplicates " + " I="+ tmp_sort_bam + " O="+ outfile + " M="+ tmp_markdup_mat + " BARCODE_TAG="+"BC" + " 2>/dev/null"
	run_cmd(cmd, "mark_dup_bam")

	cmd="samtools " + " index " + outfile + " -@ " + str(threads)
	run_cmd(cmd, "index_bam")

	#QC
	samtools_path = shutil.which("samtools")
	cmd="perl " + python_path + "AlignStat_WGS.pl " + " -in " + outfile + " -sam " + filename + " -samtools " + samtools_path + " -outdir " + outdir + " -move " + outdir
	run_cmd(cmd, "stat_bam")

	#remove tmp files
	os.remove(tmp_b_fq1)
	os.remove(tmp_b_fq2)
	os.remove(tmp_u_fq1)
	os.remove(tmp_u_fq2)
	os.remove(tmp_ema_bam)
	os.remove(tmp_bwa_bam)
	os.remove(tmp_merge_bam)
	os.remove(tmp_sort_bam)


def align_stLFR(fq1, fq2, rg, ref, outfile, sort, mark, threads):
	###########################################################################
	# To align clean ULRF (bc=30bp) sequencing reads to reference genome      #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir,filename)=os.path.split(outdir + "/" + outfile)

	#temp paths and files
	tmp_b_fq1 = outdir + "/tmp.barcoded.1.fq"
	tmp_b_fq2 = outdir + "/tmp.barcoded.2.fq"
	tmp_u_fq1 = outdir + "/tmp.unbarcoded.1.fq"
	tmp_u_fq2 = outdir + "/tmp.unbarcoded.2.fq"
	tmp_ema_sam  = outdir + "/tmp.ema.sam"
	tmp_ema_bam  = outdir + "/tmp.ema.bam"
	tmp_bwa_sam  = outdir + "/tmp.bwa.sam"
	tmp_bwa_bam  = outdir + "/tmp.bwa.bam"
	tmp_merge_bam = outdir + "/tmp.merge.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"
	rg="\'"+rg+"\'"

	#split and sort reads
	cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq1 + " | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_b_fq1
	cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq2 + " | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_b_fq2
	cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq1 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_u_fq1
	cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq2 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_u_fq2
	run_cmd(cmd1, "sort_barcoded_FQ1")
	run_cmd(cmd2, "sort_unbarcoded_FQ1")
	run_cmd(cmd3, "sort_barcoded_FQ2")
	run_cmd(cmd4, "sort_unbarcoded_FQ2")

	#Align reads
	cmd="chmod 755 " + python_path + "EMA/ema"
	subprocess.call(cmd, shell=True)
	cmd=python_path + "EMA/ema"+" align "+" -1 "+ tmp_b_fq1 + " -2 "+ tmp_b_fq2 +" -r "+ ref + " -R "+ rg + " -p "+ " stlfr "+" -t "+ str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_ema_bam + " 2>/dev/null"
	run_cmd(cmd, "EMA_alignment")
	cmd="bwa"+" mem "+" -1 " + ref + " " + tmp_u_fq1 + " "+ tmp_u_fq2 + " -R "+ rg + " -t "+ str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_bwa_bam + " 2>/dev/null"
	run_cmd(cmd, "BWA_alignment")

	#Bam Post Processing
	cmd="samtools "+ "merge " + tmp_merge_bam + " " + tmp_bwa_bam + " " + tmp_ema_bam + " --threads " + str(threads)
	run_cmd(cmd, "merge_bam")
	cmd="samtools "+ "sort " + tmp_merge_bam + " -o " + tmp_sort_bam + " --threads " + str(threads)
	run_cmd(cmd, "sort_bam")
	cmd="picard " + " MarkDuplicates " + " I= " + tmp_sort_bam + " O= " + outfile + " M= " + tmp_markdup_mat + " BARCODE_TAG= " + " BC 2>/dev/null"
	run_cmd(cmd, "Mark_bam")
	cmd="samtools " + " index " + outfile + " -@ " + str(threads)
	run_cmd(cmd, "index_bam")

	#QC
	samtools_path = shutil.which("samtools")
	cmd="perl " + python_path + "AlignStat_WGS.pl " + "-in " + outfile +" -sam " + filename + " -samtools " + samtools_path + " -outdir " + outdir + " -move " + outdir
	run_cmd(cmd,"bam_stat")

	#remove tmp files
	os.remove(tmp_b_fq1)
	os.remove(tmp_b_fq2)
	os.remove(tmp_u_fq1)
	os.remove(tmp_u_fq2)
	os.remove(tmp_bwa_bam)
	os.remove(tmp_ema_bam)
	os.remove(tmp_merge_bam)
	os.remove(tmp_sort_bam)

def align_TELLSeq(fq1, fq2, rg, ref, outfile, sort, mark, threads):
	###########################################################################
	# To align clean ULRF (bc=18bp) sequencing reads to reference genome      #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir,filename)=os.path.split(outdir + "/" + outfile)

	#temp paths and files
	tmp_b_fq1 = outdir + "/tmp.barcoded.1.fq"
	tmp_b_fq2 = outdir + "/tmp.barcoded.2.fq"
	tmp_u_fq1 = outdir + "/tmp.unbarcoded.1.fq"
	tmp_u_fq2 = outdir + "/tmp.unbarcoded.2.fq"
	tmp_ema_sam   = outdir + "/tmp.ema.sam"
	tmp_ema_bam   = outdir + "/tmp.ema.bam"
	tmp_bwa_sam   = outdir + "/tmp.bwa.sam"
	tmp_bwa_bam   = outdir + "/tmp.bwa.bam"
	tmp_merge_bam = outdir + "/tmp.merge.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"
	rg="'"+rg+"'"

	#split and sort reads
	cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq1 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_b_fq1
	cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq2 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_b_fq2
	cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq1 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_u_fq1
	cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + fq2 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_u_fq2
	run_cmd(cmd1, "sort_barcoded_FQ1")
	run_cmd(cmd2, "sort_unbarcoded_FQ1")
	run_cmd(cmd3, "sort_barcoded_FQ2")
	run_cmd(cmd4, "sort_unbarcoded_FQ2")

	#Align reads
	cmd = "chmod 755 " + python_path + "EMA/ema"
	subprocess.call(cmd, shell=True)
	cmd=python_path + "EMA/ema"+" align " + " -1 " + tmp_b_fq1 + " -2 " + tmp_b_fq2 + " -r " + ref + " -R " + rg + " -p "+ " tellseq " + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_ema_bam + " 2>/dev/null"
	run_cmd(cmd, "ema_alignment")
	cmd="bwa" + " mem " + ref + " " + tmp_u_fq1 + " " + tmp_u_fq2 + " -R "+ rg + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_bwa_bam + " 2>/dev/null"
	run_cmd(cmd, "bwa_alignment")

	#Bam Post Processing
	cmd="samtools " + " merge " + tmp_merge_bam + " " + tmp_bwa_bam + " " + tmp_ema_bam + " --threads " + str(threads) 
	run_cmd(cmd, "merge_bam")

	cmd="samtools" + " sort " + tmp_merge_bam + " -o " + tmp_sort_bam + " --threads " + str(threads) 
	run_cmd(cmd, "sort_bam")

	cmd ="picard " + " MarkDuplicates " + " I= " + tmp_sort_bam + " O= " + outfile + " M= " + tmp_markdup_mat + " BARCODE_TAG=" + "BC" + " 2>/dev/null"
	run_cmd(cmd, "mark_dup_bam")

	cmd="samtools " + " index " + outfile + " -@ " + str(threads)
	run_cmd(cmd, "index_bam")

	#QC
	samtools_path=shutil.which("samtools")
	cmd="perl " + python_path + "AlignStat_WGS.pl " + "-in " + outfile +" -sam " + filename + " -samtools " + samtools_path + " -outdir " + outdir + " -move " + outdir
	run_cmd(cmd,"bam_stat")

	#remove tmp files
	os.remove(tmp_b_fq1)
	os.remove(tmp_b_fq2)
	os.remove(tmp_u_fq1)
	os.remove(tmp_u_fq2)
	os.remove(tmp_bwa_bam)
	os.remove(tmp_ema_bam)
	os.remove(tmp_merge_bam)
	os.remove(tmp_sort_bam)

def metagenome_align(fq1, fq2, rg, ref, outfile, sort, mark, platform, threads):
	###########################################################################
	# To align clean linked read sequencing reads to reference meta genome    #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir,filename)=os.path.split(outdir + "/" + outfile)

	if (platform == "stLFR"):
		platform = "stlfr"
	elif (platform == "TELLSeq"): 
		platform = "tellseq"

	#temp paths and files
	tmp_b_fq1 = outdir + "/tmp.barcoded.1.fq"
	tmp_b_fq2 = outdir + "/tmp.barcoded.2.fq"
	tmp_u_fq1 = outdir + "/tmp.unbarcoded.1.fq"
	tmp_u_fq2 = outdir + "/tmp.unbarcoded.2.fq"
	tmp_ema_bam  = outdir + "/tmp.ema.bam"
	tmp_bwa_bam  = outdir + "/tmp.bwa.bam"
	tmp_merge_bam = outdir + "/tmp.merge.bam"
	tmp_sort_bam  = outdir + "/tmp.sort.bam"
	tmp_markdup_bam = outdir + "/tmp.markdup.bam"
	tmp_markdup_mat = outdir + "/tmp.markdup.mat"
	bam1 = outdir + "/align.round1.bam"
	bam2 = outfile
	reftxt = ref.replace("fa","txt",1)
	rg="'"+rg+"'"

	#split and sort reads
	cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"|\")}\' " + fq1 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \" \" |tr \"|\" \"\\n\"  > " + tmp_b_fq1
	cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"|\")}\' " + fq2 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \" \" |tr \"|\" \"\\n\" > " + tmp_b_fq2
	cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"|\")}\' " + fq1 + " | grep -v \"BX:Z:\" | tr \"|\" \"\\n\" > " + tmp_u_fq1
	cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"|\")}\' " + fq2 + " | grep -v \"BX:Z:\" | tr \"|\" \"\\n\" > " + tmp_u_fq2
	run_cmd(cmd1, "sort_barcode_FQ1")
	run_cmd(cmd2, "sort_barcode_FQ2")
	run_cmd(cmd3, "sort_unbarcode_FQ1")
	run_cmd(cmd4, "sort_unbarcode_FQ2")

	#Align reads round1
	cmd = "chmod 755 " + python_path + "EMA/ema"
	subprocess.call(cmd, shell=True)
	cmd=python_path + "EMA/ema"+" align "+" -1 "+ tmp_b_fq1 + " -2 "+ tmp_b_fq2 +" -r "+ ref + " -R "+ rg + " -p "+ platform + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_ema_bam + " --threads " + str(threads) + " 2>/dev/null"
	run_cmd(cmd, "ema_align")
	cmd="bwa mem " + ref + " " + tmp_u_fq1 + " " + tmp_u_fq2 + " -R "+ rg + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " -o " + tmp_bwa_bam  + " --threads " + str(threads) + " 2>/dev/null"
	run_cmd(cmd, "bwa_align")
	cmd="samtools merge " + tmp_merge_bam + " " + tmp_bwa_bam + " " + tmp_ema_bam + " --threads " + str(threads)
	subprocess.call(cmd, shell=True)
	cmd="samtools sort " + " --threads " + str(threads) + " " + tmp_merge_bam + " -o " + bam1
	subprocess.call(cmd, shell=True)
	cmd="samtools index " + bam1 + " -@  " + str(threads)
	subprocess.call(cmd, shell=True)
	cmd="bedtools makewindows -g " + reftxt + " -w 1000 " + " > " + outdir + "/windows.bed"
	run_cmd(cmd, "makewindows")
	cmd="perl " + python_path + "Readcount.Bam2bed.pl" + " -db " + outdir + "/windows.bed" + " -bam " + bam1 + " -outfile " + outdir + "/windows.depth.xls"
	run_cmd(cmd, "readcount")
	cmd="python " + python_path + "CalculateCoverage.round1.py " + outdir + "/windows.depth.xls " + outdir + "/UHGG.coverage_average.xls"
	run_cmd(cmd, "cal_cov")
	cmd="rm -rf " + outdir + "/tmp.merge.bam " + outdir + "/tmp.bwa.bam " + outdir + "/tmp.ema.bam" + outdir + "/align.round1.bam" + outdir + "/align.round1.bam.bai"
	subprocess.call(cmd, shell=True)

	#Align reads round2
	cmd="python " + python_path + "ExtractRef.py " + outdir + "/UHGG.coverage_average.xls " + outdir + "/UHGG.reference.subset.fa " + outdir + "/UHGG.reference.subset.txt"
	subprocess.call(cmd, shell=True)
	genome2=outdir + "/UHGG.reference.subset.txt"
	cmd="bwa index " + outdir + "/UHGG.reference.subset.fa" + " 2>/dev/null"
	subprocess.call(cmd, shell=True)
	cmd="samtools faidx " + outdir + "/UHGG.reference.subset.fa" + " 2>/dev/null"
	subprocess.call(cmd, shell=True)
	cmd=python_path + "EMA/ema"+" align "+" -1 "+ tmp_b_fq1 + " -2 "+ tmp_b_fq2 +" -r "+ outdir + "/UHGG.reference.subset.fa" + " -R " + rg + " -p "+ platform + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " - " + " -o " + tmp_ema_bam + " --threads " + str(threads) + " 2>/dev/null"
	run_cmd(cmd, "ema_align")
	cmd="bwa mem " + outdir + "/UHGG.reference.subset.fa" + " " + tmp_u_fq1 + " " + tmp_u_fq2 + " -R "+ rg + " -t " + str(threads) + " | " + "samtools" + " view " + " -Sb " + " - " + " -o " + tmp_bwa_bam  + " --threads " + str(threads) + " 2>/dev/null"
	run_cmd(cmd, "bwa_align")
	cmd="samtools merge -f " + tmp_merge_bam + " " + tmp_bwa_bam + " " + tmp_ema_bam + " --threads " + str(threads)
	subprocess.call(cmd, shell=True)
	cmd="samtools sort " + " --threads " + str(threads) + " " + tmp_merge_bam + " -o " + bam2
	subprocess.call(cmd, shell=True)
	cmd="samtools index " + bam2 + " -@  " + str(threads)
	subprocess.call(cmd, shell=True)
	cmd="bedtools makewindows -g " + outdir + "/UHGG.reference.subset.txt" + " -w 1000 " + " > " + outdir + "/windows.iter.bed"
	run_cmd(cmd, "makewindows")
	cmd="perl " + python_path + "Readcount.Bam2bed.pl" + " -db " + outdir + "/windows.iter.bed" + " -bam " + bam2 + " -outfile " + outdir + "/windows.depth.iter.xls"
	run_cmd(cmd, "readcount")
	cmd="python " + python_path + "CalculateCoverage.round2.py " + outdir + "/windows.depth.iter.xls " + outdir + "/UHGG.coverage_average.iter.xls"
	run_cmd(cmd, "cal_cov")
	cmd="rm -rf " + outdir + "/tmp.merge.bam " + outdir + "/tmp.bwa.bam " + outdir + "/tmp.ema.bam"
	subprocess.call(cmd, shell=True)

	#remove tmp files
	os.remove(tmp_b_fq1)
	os.remove(tmp_b_fq2)
	os.remove(tmp_u_fq1)
	os.remove(tmp_u_fq2)
	os.remove(outdir + "/windows.iter.bed")
	os.remove(outdir + "/windows.depth.iter.xls")
	os.remove(outdir + "/windows.bed")
	os.remove(outdir + "/windows.depth.xls")

