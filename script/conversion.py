import os
import subprocess
import sys

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

from utility import *


def TENx2ULRF(infq1, infq2, outfq1, outfq2, white_list, host, filt, sort, threads, genome):
	###########################################################################
	# To convert the raw 10x sequencing data to unified linked read format.   #
	###########################################################################	
	if outfq1.startswith("/"):
		(outdir, filename1) = os.path.split(outfq1)
		(outdir, filename2) = os.path.split(outfq2)
	else:
		outdir = os.getcwd()
		filename1 = outfq1
		filename2 = outfq2

	#temp paths and files
	bc_fq  = outdir + "/tmp_bc.fq"
	bc_sai = outdir + "/tmp_bc.sai"
	bc_sam = outdir + "/tmp_bc.sam"
	tmp_fq1_read = outdir + "/tmp.read1_only.fq"
	tmp_out1 = outdir + "/tmp.read1_out.fq"
	tmp_out2 = outdir + "/tmp.read2_out.fq"
	tmp_filt1 = outdir + "/tmp.read1.filt.fq"
	tmp_filt2 = outdir + "/tmp.read2.filt.fq"
	tmp_filt_sort_wb1  = outdir + "/" + filename1 + ".sort.wb.fq"
	tmp_filt_sort_wob1 = outdir + "/" + filename1 + ".sort.wob.fq"
	tmp_filt_sort_wb2  = outdir + "/" + filename2 + ".sort.wb.fq"
	tmp_filt_sort_wob2 = outdir + "/" + filename2 + ".sort.wob.fq"

	reporthtml = outdir + "/" + "FASTQ.QC.html"
	reportjson = outdir + "/" + "FASTQ.QC.json"

	#extract barcodes
	cmd="python "+python_path+"extract_bc_10x.py "+"--infile "+infq1+" --outfile "+tmp_fq1_read+" --outbarcode "+bc_fq
	run_cmd(cmd,"extract_bc")
	#align barcodes
	cmd="bwa aln " + "-t " + str(threads) + " " + white_list + " " + bc_fq + " -f " + bc_sai + " 2>/dev/null"
	run_cmd(cmd,"aln_bc")
	cmd="bwa "+ "samse "+ white_list + " " + bc_sai + " " + bc_fq + " -f " + bc_sam + " 2>/dev/null"
	run_cmd(cmd, "samse_bc")
	#correct barcodes
	cmd="python " + python_path + "correct_bc_10x.py " + "-1 " + tmp_fq1_read + " -2 " + infq2 + " -w " + white_list + " -s " + bc_sam + " -3 " + tmp_out1 + " -4 " + tmp_out2 + " -n " + str("2")
	run_cmd(cmd, "correct_bc_10x")

	if ((sort == "Yes") & (filt == "Yes")):
		#filtering reads start
		cmd="fastp "+"-i "+tmp_out1+" -I "+tmp_out2+" -o "+tmp_filt1+" -O "+tmp_filt2+" -w 8 -q 20 -u 50 -n 10 -l 50 -h "+reporthtml+" -j " + reportjson + " 2>/dev/null"
		run_cmd(cmd, "Filt_FQ")

		#sorting reads by barcodes start
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		run_cmd(cmd1, "sort_barcoded_FQ1")
		run_cmd(cmd2, "sort_unbarcoded_FQ1")

		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		run_cmd(cmd3, "sort_barcoded_FQ2")
		run_cmd(cmd4, "sort_unbarcoded_FQ2")

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
		#logging("sorting reads by barcodes end.")

		os.remove(tmp_out1)
		os.remove(tmp_out2)
		os.remove(tmp_filt1)
		os.remove(tmp_filt2)
		os.remove(tmp_filt_sort_wb1)
		os.remove(tmp_filt_sort_wb2)
		os.remove(tmp_filt_sort_wob1)
		os.remove(tmp_filt_sort_wob2)
	elif (filt == "Yes"):
		#filtering reads start
		cmd="fastp "+"-i "+tmp_out1+" -I "+tmp_out2+" -o "+outfq1+" -O "+outfq2+" -w 8 -q 20 -u 50 -n 10 -l 50 -h "+reporthtml+" -j "+reportjson + " 2>/dev/null" 
		run_cmd(cmd, "Filt_FQ")
		os.remove(tmp_out1)
		os.remove(tmp_out2)
	elif (sort == "Yes"):
		#sorting reads by barcodes start
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		run_cmd(cmd1, "sort_barcoded_FQ1")
		run_cmd(cmd2, "sort_unbarcoded_FQ1")

		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		run_cmd(cmd3, "sort_barcoded_FQ2")
		run_cmd(cmd4, "sort_unbarcoded_FQ2")

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)

		os.remove(tmp_filt_sort_wb1)
		os.remove(tmp_filt_sort_wb2)
		os.remove(tmp_filt_sort_wob1)
		os.remove(tmp_filt_sort_wob2)
	else:
		subprocess.call(["mv", tmp_out1, outfq1])
		subprocess.call(["mv", tmp_out2, outfq2])

	#metagenomic specific section
	if(genome == "metagenome"):
		rmhost_fq1=outdir + "/tmp.rmhost.1.fq"
		rmhost_fq2=outdir + "/tmp.rmhost.2.fq"
		hostbam=outdir + "/tmp.HUMAN.align.bam"
		hostunmapped=outdir + "/tmp.HUMAN.unmapped.bam"

		cmd1="bwa mem -C " + host + " " + outfq1 + " " + outfq2 + " -t " + str(threads) + " |samtools sort -@ " + str(threads) + " -o " + hostbam + " 2>/dev/null"
		run_cmd(cmd1, "Host_alignment")
		cmd2="samtools view -b -@ " + str(threads) + " -f 4 " + hostbam + " -o " + hostunmapped
		run_cmd(cmd2, "Extract_unmapped_read")
		cmd3="samtools fastq -@ " + str(threads) + " " + hostunmapped + " -1 " + rmhost_fq1 + " -2 " + rmhost_fq2 + " -T BX -s /dev/null"
		run_cmd(cmd3, "BAM2FQ")
		subprocess.call(["mv", rmhost_fq1, outfq1])
		subprocess.call(["mv", rmhost_fq2, outfq2])
		os.remove(hostbam)
		os.remove(hostunmapped)

	#remove tmp files
	os.remove(tmp_fq1_read)
	os.remove(bc_sam)
	os.remove(bc_sai)
	os.remove(bc_fq)


def stLFR2ULRF(infq1, infq2, outfq1, outfq2, white_list, host, filt, sort, threads, genome):
	###########################################################################
	# To convert the raw stLFR sequencing data to unified linked read format  #
	###########################################################################

	if outfq1.startswith("/"):
		(outdir, filename1) = os.path.split(outfq1)
		(outdir, filename2) = os.path.split(outfq2)
	else:
		outdir = os.getcwd()
		filename1=outfq1
		filename2=outfq2

	#temp paths and files
	bc_fq   = outdir + "/tmp"
	bc1_fq  = outdir + "/tmp_bc1.fq"
	bc2_fq  = outdir + "/tmp_bc2.fq"
	bc3_fq  = outdir + "/tmp_bc3.fq"

	bc1_sai = outdir + "/tmp_bc1.sai"
	bc2_sai = outdir + "/tmp_bc2.sai"
	bc3_sai = outdir + "/tmp_bc3.sai"

	bc1_sam = outdir + "/tmp_bc1.sam"
	bc2_sam = outdir + "/tmp_bc2.sam"
	bc3_sam = outdir + "/tmp_bc3.sam"
	
	tmp_fq2_read = outdir + "/tmp_read2_only.fq"
	tmp_out1 = outdir + "/tmp_read1_out.fq"
	tmp_out2 = outdir + "/tmp_read2_out.fq"
	tmp_filt1 = outdir + "/tmp.read1.filt.fq"
	tmp_filt2 = outdir + "/tmp.read2.filt.fq"
	tmp_filt_sort_wb1  = outdir + "/" + filename1 + ".sort.wb.fq"
	tmp_filt_sort_wob1 = outdir + "/" + filename1 + ".sort.wob.fq"
	tmp_filt_sort_wb2  = outdir + "/" + filename2 + ".sort.wb.fq"
	tmp_filt_sort_wob2 = outdir + "/" + filename2 + ".sort.wob.fq"
	
	reporthtml = outdir + "/" + "FASTQ.QC.html"
	reportjson = outdir + "/" + "FASTQ.QC.json"

	#extract barcodes
	cmd="python "+python_path+"extract_bc_stLFR.py "+infq2+" "+bc_fq
	run_cmd(cmd,"extract_bc")

	#align barcodes
	cmd="bwa aln -t "+str(threads)+" "+white_list+" "+bc1_fq+" -f "+ bc1_sai + " 2>/dev/null"
	run_cmd(cmd,"aln_bc1")
	cmd="bwa samse "+white_list+" "+bc1_sai+" "+bc1_fq+" -f "+ bc1_sam + " 2>/dev/null"
	run_cmd(cmd,"samse_bc1")
	cmd="bwa aln -t "+str(threads)+" "+white_list+" "+bc2_fq+" -f "+ bc2_sai + " 2>/dev/null"
	run_cmd(cmd, "aln_bc2")
	cmd="bwa samse "+white_list+" "+bc2_sai+" "+bc2_fq+" -f "+bc2_sam + " 2>/dev/null"
	run_cmd(cmd, "samse_bc2")
	cmd="bwa aln -t "+str(threads)+" "+white_list+" "+bc3_fq+" -f "+bc3_sai + " 2>/dev/null"
	run_cmd(cmd, "aln_bc3")
	cmd="bwa samse "+white_list+" "+bc3_sai+" "+bc3_fq+" -f "+bc3_sam + " 2>/dev/null"
	run_cmd(cmd, "samse_bc3")

	#correct barcodes
	cmd="chmod 755 " + python_path + "/correct_barcode_stlfr"
	subprocess.call(cmd, shell=True)
	cmd=python_path+"/correct_barcode_stlfr "+white_list+" "+bc1_sam+" "+bc2_sam+" "+bc3_sam+" "+infq1+" "+tmp_fq2_read+" "+str(2)+" "+tmp_out1+" "+tmp_out2
	run_cmd(cmd, "barcode_correction")

	#Filt & sort barcodes
	if ((sort == "Yes") & (filt == "Yes")):
		cmd="fastp -i "+tmp_out1+" -I "+tmp_out2+" -o "+tmp_filt1+" -O "+tmp_filt2+" -w 8 -q 20 -u 50 -n 10 -l 50 -h "+reporthtml+" -j "+reportjson+ " 2>/dev/null" 
		run_cmd(cmd, "filt_FQ")

		#sorting reads by barcodes start
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + "  | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		run_cmd(cmd1, "sort_barcoded_FQ1")
		run_cmd(cmd2, "sort_unbarcoded_FQ1")

		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + "  | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		run_cmd(cmd3, "sort_barcoded_FQ2")
		run_cmd(cmd4, "sort_unbarcoded_FQ2")

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
		os.remove(tmp_out1)
		os.remove(tmp_out2)
		os.remove(tmp_filt1)
		os.remove(tmp_filt2)
		os.remove(tmp_filt_sort_wb1)
		os.remove(tmp_filt_sort_wb2)
		os.remove(tmp_filt_sort_wob1)
		os.remove(tmp_filt_sort_wob2)
		#logging("sorting reads by barcodes end.")
	elif (filt == "Yes"):
		cmd="fastp "+"-i "+tmp_out1+" -I "+tmp_out2+" -o "+outfq1+" -O "+outfq2+" -w 8 -q 20 -u 50 -n 10 -l 50 -h "+reporthtml+" -j "+reportjson + " 2>/dev/null" 
		run_cmd(cmd, "filt_FQ")

		os.remove(tmp_out1)
		os.remove(tmp_out2)
	elif (sort == "Yes"):
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		run_cmd(cmd1, "sort_barcoded_FQ1")
		run_cmd(cmd2, "sort_unbarcoded_FQ1")
		
		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		run_cmd(cmd3, "sort_barcoded_FQ2")
		run_cmd(cmd4, "sort_unbarcoded_FQ2")

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)

		os.remove(tmp_out1)
		os.remove(tmp_out2)
		os.remove(tmp_filt_sort_wb1)
		os.remove(tmp_filt_sort_wb2)
		os.remove(tmp_filt_sort_wob1)
		os.remove(tmp_filt_sort_wob2)
	else:
		subprocess.call(["mv", tmp_out1, outfq1])
		subprocess.call(["mv", tmp_out2, outfq2])

	#metagenomic specific section
	if(genome == "metagenome"):
		rmhost_fq1=outdir + "/tmp.rmhost.1.fq"
		rmhost_fq2=outdir + "/tmp.rmhost.2.fq"
		hostbam=outdir + "/tmp.HUMAN.align.bam"
		hostunmapped=outdir + "/tmp.HUMAN.unmapped.bam"

		cmd1="bwa mem -C " + host + " " + outfq1 + " " + outfq2 + " -t " + str(threads) + " |samtools sort -@ " + str(threads) + " -o " + hostbam + " 2>/dev/null"
		run_cmd(cmd1, "Host_alignment")
		cmd2="samtools view -b -@ "+str(threads)+" -f 4 "+hostbam+" -o "+hostunmapped
		run_cmd(cmd2, "Extract_unmapped_read")
		cmd3="samtools fastq -@ "+str(threads)+" "+hostunmapped+" -1 "+rmhost_fq1+" -2 "+rmhost_fq2+"  -T BX -s /dev/null"
		run_cmd(cmd3, "BAM2FQ")
		subprocess.call(["mv", rmhost_fq1, outfq1])
		subprocess.call(["mv", rmhost_fq2, outfq2])
		os.remove(hostbam)
		os.remove(hostunmapped)

	#remove tmp files
	os.remove(tmp_fq2_read)
	os.remove(bc1_fq)
	os.remove(bc2_fq)
	os.remove(bc3_fq)
	os.remove(bc1_sai)
	os.remove(bc2_sai)
	os.remove(bc3_sai)
	os.remove(bc1_sam)
	os.remove(bc2_sam)
	os.remove(bc3_sam)

def TELLSeq2ULRF(infq1, infq2, idxfq, outfq1, outfq2, host, filt, sort, threads, genome):
	###########################################################################
	# To convert raw TELL-Seq sequencing data to univfied linked read format. #
	###########################################################################
	if outfq1.startswith("/"):
		(outdir, filename1) = os.path.split(outfq1)
		(outdir, filename2) = os.path.split(outfq2)
	else:
		outdir = os.getcwd()
		filename1=outfq1
		filename2=outfq2

	#temp files
	tmpbq   = outdir + "/tmp.barcode.txt"
	tmpcbq  = outdir + "/tmp.barcode.correct.txt"
	tmpebq  = outdir + "/tmp.barcode.error.txt"
	tmpfq1  = outdir + "/tmp.correct.1.fq"
	tmpfq2  = outdir + "/tmp.correct.2.fq"	
	tmp_filt1 = outdir + "/tmp.read1.filt.fq"
	tmp_filt2 = outdir + "/tmp.read2.filt.fq"
	tmp_filt_sort_wb1  = outdir + "/" + filename1 + ".sort.wb.fq"
	tmp_filt_sort_wob1 = outdir + "/" + filename1 + ".sort.wob.fq"
	tmp_filt_sort_wb2  = outdir + "/" + filename2 + ".sort.wb.fq"
	tmp_filt_sort_wob2 = outdir + "/" + filename2 + ".sort.wob.fq"
	reporthtml = outdir + "/" + "FASTQ.QC.html"
	reportjson = outdir + "/" + "FASTQ.QC.json"

	#Format conversion
	cmd="python " + python_path + "extract_bc_TELLSeq.py " + idxfq + " " + tmpbq
	run_cmd(cmd, "Extract_bc")
	cmd="python " + python_path + "filt_bc_TELLSeq.py -I " + tmpbq + " -O1 " + tmpcbq + " -O2 " + tmpebq
	run_cmd(cmd, "filt_bc")
	#correct_barcode
	cmd="python "+python_path+"correct_bc.TELLSeq.py -F1 " + infq1 + " -F2 " + infq2 + " -B " + idxfq + " -O1 " + tmpfq1 + " -O2 " + tmpfq2 + " -C " + tmpcbq
	run_cmd(cmd, "correct_barcode")

	#Filt & sort barcodes
	if ((sort == "Yes") & (filt == "Yes")):
		cmd="fastp -i "+tmpfq1+" -I "+tmpfq2+" -o "+tmp_filt1+" -O "+tmp_filt2+" -w 8 -q 20 -u 50 -n 10 -l 50 -h "+reporthtml+" -j "+reportjson + " 2>/dev/null" 
		run_cmd(cmd, "filt_FQ")

		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wb2
		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		run_cmd(cmd1, "sort_barcoded_FQ1")
		run_cmd(cmd2, "sort_barcoded_FQ2")
		run_cmd(cmd3, "sort_unbarcoded_FQ1")
		run_cmd(cmd4, "sort_unbarcoded_FQ2")

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)

		os.remove(tmpfq1)
		os.remove(tmpfq2)
		os.remove(tmp_filt1)
		os.remove(tmp_filt2)
		os.remove(tmp_filt_sort_wb1)
		os.remove(tmp_filt_sort_wb2)
		os.remove(tmp_filt_sort_wob1)
		os.remove(tmp_filt_sort_wob2)

	elif (filt == "Yes"):
		cmd="fastp "+"-i "+tmpfq1+" -I "+tmpfq2+" -o "+outfq1+" -O "+outfq2+" -w 8 -q 20 -u 50 -n 10 -l 50 -h "+reporthtml+" -j "+reportjson + " 2>/dev/null"
		run_cmd(cmd, "filt_FQ")

		os.remove(tmpfq1)
		os.remove(tmpfq2)
	elif (sort == "Yes"):
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wob1
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep -v \"BX:Z:\" | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wob2
		run_cmd(cmd1, "sort_barcoded_FQ1")
		run_cmd(cmd2, "sort_unbarcoded_FQ1")
		run_cmd(cmd3, "sort_barcoded_FQ2")
		run_cmd(cmd4, "sort_unbarcoded_FQ2")

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
        
		os.remove(tmpfq1)
		os.remove(tmpfq2)
		os.remove(tmp_filt_sort_wb1)
		os.remove(tmp_filt_sort_wb2)
		os.remove(tmp_filt_sort_wob1)
		os.remove(tmp_filt_sort_wob2)
	else:
		subprocess.call(["mv", tmpfq1, outfq1])
		subprocess.call(["mv", tmpfq2, outfq2])

	#metagenomic specific section
	if(genome == "metagenome"):
		rmhost_fq1=outdir + "/tmp.rmhost.1.fq"
		rmhost_fq2=outdir + "/tmp.rmhost.2.fq"
		hostbam=outdir + "/tmp.HUMAN.align.bam"
		hostunmapped=outdir + "/tmp.HUMAN.unmapped.bam"

		cmd1="bwa mem -C " + host + " " + outfq1 + " " + outfq2 + " -t " + str(threads) + " |samtools sort -@ " + str(threads) + " -o " + hostbam + " 2>/dev/null"
		run_cmd(cmd1, "Host_alignment")
		cmd2="samtools view -b -@ "+ str(threads) + " -f 4 " + hostbam + " -o " + hostunmapped
		run_cmd(cmd2, "extract_unmapped_reads")
		cmd3="samtools fastq -@ " + str(threads) + " " + hostunmapped + " -1 " + rmhost_fq1 + " -2 " + rmhost_fq2 + " -T BX -s /dev/null"
		run_cmd(cmd3, "BAM2FQ")
		subprocess.call(["mv", rmhost_fq1, outfq1])
		subprocess.call(["mv", rmhost_fq2, outfq2])
		os.remove(hostbam)
		os.remove(hostunmapped)

	#remove tmp files
	os.remove(tmpbq)
	os.remove(tmpcbq)
	os.remove(tmpebq)
