import os
import subprocess
import sys

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

from utility import *


def TENx2ULRF(infq1, infq2, outfq1, outfq2, white_list, filt, sort, threads, bin):
	###########################################################################
	# To convert the raw 10x sequencing data to unified linked read format.   #
	###########################################################################	
	if outfq1.startswith("/"):
		(outdir, filename1) = os.path.split(outfq1)
		(outdir, filename2) = os.path.split(outfq2)
	else:
		outdir = os.getcwd()

	#temp paths and files
	bc_fq  = outdir + "/tmp_bc.fq"
	bc_sai = outdir + "/tmp_bc.sai"
	bc_sam = outdir + "/tmp_bc.sam"
	tmp_fq1_read = outdir + "/tmp.read1_only.fq"
	tmp_out1 = outdir + "/tmp.read1_out.fq"
	tmp_out2 = outdir + "/tmp.read2_out.fq"
	tmp_filt1 = outdir + "/tmp.read1.filt.fq"
	tmp_filt2 = outdir + "/tmp.read2.filt.fq"
	tmp_filt_sort_wb1  = outdir + "/" + filename1 + ".sort.wb.1.fq"
	tmp_filt_sort_wob1 = outdir + "/" + filename1 + ".sort.wob.1.fq"
	tmp_filt_sort_wb2  = outdir + "/" + filename2 + ".sort.wb.2.fq"
	tmp_filt_sort_wob2 = outdir + "/" + filename2 + ".sort.wob.2.fq"
	
	reporthtml = outdir + "/FASTQ.QC.fastp.html"
	reportjson = outdir + "/FASTQ.QC.fastp.json"

	#extract barcodes
	logging("Extracting barcodes start.")
	run_cmd([
		"python", 
		bin+"/python/extract_bc_10x.py",
		"--infile", infq1,
		"--outfile", tmp_fq1_read, 
		"--outbarcode", bc_fq
		],
		"Extract_barcodes"
		)
	logging("Extracting barcodes end.")

	#align barcodes
	logging("aligning barcodes start.")
	run_cmd(
		[
		"bwa",
	 	"aln", 
		"-t", threads, 
		white_list,
		bc_fq, 
		"-f", bc_sai
		],
		"aln_barcodes"
	)

	run_cmd(
		[
		"bwa",
		"samse", 
		white_list, 
		bc_sai, 
		bc_fq, 
		"-f", 
		bc_sam
		],
		"samse_barcodes"
	)
	logging("aligning barcodes end.")

	#correct barcodes
	logging("correcting barcodes start.")
	run_cmd(
		[
		"python", 
		bin+"/python/correct_bc_10x.py", 
		"-1", tmp_fq1_read, 
		"-2", infq2, 
		"-w", white_list, 
		"-s", bc_sam, 
		"-3", tmp_out1, 
		"-4", tmp_out2, 
		"-n", "2"
		],
		"correct_barcodes"
	)
	logging("correcting barcodes end.")

	#Filt & sort barcodes
	logging("fitering fastq start.")
	if ((sort == "Yes") & (filt == "Yes")):
		#filtering reads start
		run_cmd(
			[
			"fastp",
			"-i", tmp_out1,
			"-I", tmp_out2, 
			"-o", tmp_filt1, 
			"-O", tmp_filt2, 
			"-w", threads, 
			"-q", "20", 
			"-u","50", 
			"-n", "10", 
			"-l", "50",
			"-h", reporthtml, 
			"-j", reportjson
			],
			"filter_reads"
		)
		#filtering reads end
	
		#sorting reads by barcodes start
		logging("sorting reads by barcodes start.")
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		subprocess.call(cmd1, shell=True)
		subprocess.call(cmd2, shell=True)

		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		subprocess.call(cmd3, shell=True)
		subprocess.call(cmd4, shell=True)

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
		logging("sorting reads by barcodes end.")
		
		os.remove(tmp_out1)
		os.remove(tmp_out2)
		os.remove(tmp_filt1)
		os.remove(tmp_filt2)
		#sorting reads by barcodes end
	elif (filt == "Yes"):
		#filtering reads start
		run_cmd(
			[
			"fastp",
			"-i", tmp_out1,
			"-I", tmp_out2,
			"-o", outfq1,
			"-O", outfq2,
			"-w", threads,
			"-q", "20",
			"-u","50",
			"-n", "10",
			"-l", "50",
			"-h", reporthtml,
			"-j", reportjson
			],
			"filter_reads"
                )
                #filtering reads end
		os.remove(tmp_out1)
		os.remove(tmp_out2)		
	elif (sort == "Yes"):
		#sorting reads by barcodes start
		logging("sorting reads by barcodes start.")
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		subprocess.call(cmd1, shell=True)
		subprocess.call(cmd2, shell=True)

		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep \"BX:Z:\" | sort -k 2.1,2.21 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		subprocess.call(cmd3, shell=True)
		subprocess.call(cmd4, shell=True)

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
		logging("sorting reads by barcodes end.")
	else:
		subprocess.call(["mv", tmp_out1, outfq1])
		subprocess.call(["mv", tmp_out2, outfq2])

	#remove tmp files
	os.remove(tmp_fq1_read)
	os.remove(bc_sam)
	os.remove(bc_sai)
	os.remove(bc_fq)


def stLFR2ULRF(infq1, infq2, outfq1, outfq2, white_list, filt, sort, threads, bin):
	###########################################################################
	# To convert the raw stLFR sequencing data to unified linked read format  #
	###########################################################################

	if outfq1.startswith("/"):
		(outdir, filename1) = os.path.split(outfq1)
		(outdir, filename2) = os.path.split(outfq2)
	else:
		outdir = os.getcwd()

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
	
	tmp_fq1_read = outdir + "/tmp_read1_only.fq"
	tmp_out1 = outdir + "/tmp_read1_out.fq"
	tmp_out2 = outdir + "/tmp_read2_out.fq"
	tmp_filt1 = outdir + "/tmp.read1.filt.fq"
	tmp_filt2 = outdir + "/tmp.read2.filt.fq"
	tmp_filt_sort_wb1  = outdir + "/" + filename1 + ".sort.wb.1.fq"
	tmp_filt_sort_wob1 = outdir + "/" + filename1 + ".sort.wob.1.fq"
	tmp_filt_sort_wb2  = outdir + "/" + filename2 + ".sort.wb.2.fq"
	tmp_filt_sort_wob2 = outdir + "/" + filename2 + ".sort.wob.2.fq"
	
	reporthtml = outdir + "/FASTQ.QC.fastp.html"
	reportjson = outdir + "/FASTQ.QC.fastp.json"

	#extract barcodes
	logging("Extracting barcodes start.")
	run_cmd(
		[
		"python",
		bin+"/python/extract_bc_stLFR.py",
		infq2, 
		bc_fq],
		"Extract_barcode"
	)
	logging("Extracting barcodes end.")

	#align barcodes
	logging("Aligning barcodes start.")
	run_cmd(
		["bwa",
		"aln", 
		"-t", threads,
		white_list, 
		bc1_fq, 
		"-f", bc1_sai
		],
		"aln_barcode1"
		)
	
	run_cmd(
		[
		"bwa", 
		"samse", 
		white_list, 
		bc1_sai, 
		bc1_fq, 
		"-f", 
		bc1_sam
		],
		"samse_barcode1"
	)

	run_cmd(
		[
		"bwa", 
		"aln", 
		"-t", threads, 
		white_list, 
		bc2_fq, 
		"-f", bc2_sai
		],
		"aln_barcode2"
	)
	
	run_cmd(
		[
		"bwa", 
		"samse",
		white_list, 
		bc2_sai, 
		bc2_fq, 
		"-f", bc2_sam
		],
		"samse_barcode2"
	)

	run_cmd(
		[
		"bwa", 
		"aln", 
		"-t", threads,
		white_list, 
		bc3_fq, 
		"-f", bc3_sai
		],
		"aln_barcode3"
	)
	
	run_cmd(
		[
		"bwa", 
		"samse", 
		white_list, 
		bc3_sai, 
		bc3_fq, 
		"-f", bc3_sam
		],
		"samse_barcode3"
	)
	logging("Extracting barcodes end.")

	#correct barcodes
	logging("Correcting barcodes start.")
	run_cmd(
		[
		python_path + "../src/correct_barcode_stlfr",
		white_list,
		bc1_sam,
		bc2_sam,
		bc3_sam,
		tmp_fq1_read,
		infq2,
		str(2),
		tmp_out1, 
		tmp_out2,
		],
		"Correcting_barcodes"
	)
	logging("Correcting barcodes end.")

	#Filt & sort barcodes
	logging("fitering fastq start.")
	if ((sort == "Yes") & (filt == "Yes")):
		run_cmd(
			[
			"fastp",
			"-i", tmp_out1, 
			"-I", tmp_out2, 
			"-o", tmp_filt1, 
			"-O", tmp_filt2, 
			"-w", threads, 
			"-q", "20", 
			"-u","50", 
			"-n", "10", 
			"-l", "50",
			"-h", reporthtml, 
			"-j", reportjson
			],
			"Fitering_reads"
		)
		#sorting reads by barcodes start
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + "  | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		subprocess.call(cmd1, shell=True)
		subprocess.call(cmd2, shell=True)

		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + "  | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		subprocess.call(cmd3, shell=True)
		subprocess.call(cmd4, shell=True)

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
		os.remove(tmp_out1)
		os.remove(tmp_out2)
		os.remove(tmp_filt1)
		os.remove(tmp_filt2)		
		logging("sorting reads by barcodes end.")
	elif (filt == "Yes"):
		run_cmd(
			[
			"fastp", 
			"-i", tmp_out1, 
			"-I", tmp_out2, 
			"-o", outfq1, 
			"-O", outfq2, 
			"-w", threads, 
			"-q", "20", 
			"-u","50", 
			"-n", "10", 
			"-l", "50",
			"-h", reporthtml, 
			"-j", reportjson
			],
			"Fitering_reads"
		)
		
		os.remove(tmp_out1)
		os.remove(tmp_out2)
	elif (sort == "Yes"):
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out1 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob1
		subprocess.call(cmd1, shell=True)
		subprocess.call(cmd2, shell=True)
		
		cmd3 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep \"BX:Z:\" | sort -k 2.1,2.36 -k 1,1 | tr \"\\t\" \"\\n\" > " + tmp_filt_sort_wb2
		cmd4 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_out2 + " | grep -v \"BX:Z:\" | sort -k 1,1 | tr \"\\t\" \"\\n\" >" + tmp_filt_sort_wob2
		subprocess.call(cmd3, shell=True)
		subprocess.call(cmd4, shell=True)

		cmd5 = "cat " + tmp_filt_sort_wb1 + " " + tmp_filt_sort_wob1 + " > " + outfq1
		cmd6 = "cat " + tmp_filt_sort_wb2 + " " + tmp_filt_sort_wob2 + " > " + outfq2
		subprocess.call(cmd5, shell=True)
		subprocess.call(cmd6, shell=True)
		os.remove(tmp_out1)
		os.remove(tmp_out2)
	else:
		subprocess.call(["mv", tmp_out1, outfq1])
		subprocess.call(["mv", tmp_out2, outfq2])
	
	#remove tmp files
	os.remove(tmp_fq1_read)
	os.remove(bc1_fq)
	os.remove(bc2_fq)
	os.remove(bc3_fq)
	os.remove(bc1_sai)
	os.remove(bc2_sai)
	os.remove(bc3_sai)
	os.remove(bc1_sam)
	os.remove(bc2_sam)
	os.remove(bc3_sam)


def TELLSeq2ULRF(infq1, infq2, idxfq, outfq1, outfq2, filt, sort, threads, bin):
	###########################################################################
	# To convert raw TELL-Seq sequencing data to univfied linked read format. #
	###########################################################################
	if outfq1.startswith("/"):
		(outdir, filename) = os.path.split(outfq1)
	else:
		outdir = os.getcwd()

	#temp files
	tmpbq   = outdir + "/tmp.barcode.txt"
	tmpcbq  = outdir + "/tmp.barcode.correct.txt"
	tmpebq  = outdir + "/tmp.barcode.error.txt"
	tmpfq1  = outdir + "/tmp.correct.1.fq"
	tmpfq2  = outdir + "/tmp.correct.2.fq"	
	tmp_filt1 = outdir + "/tmp.read1.filt.fq"
	tmp_filt2 = outdir + "/tmp.read2.filt.fq"
	tmp_filt_sort_wb1  = outdir + "/" + filename + ".sort.wb.1.fq"
	tmp_filt_sort_wob1 = outdir + "/" + filename + ".sort.wob.1.fq"
	tmp_filt_sort_wb2  = outdir + "/" + filename + ".sort.wb.2.fq"
	tmp_filt_sort_wob2 = outdir + "/" + filename + ".sort.wob.2.fq"
	reporthtml = outdir + "/FASTQ.QC.fastp.html"
	reportjson = outdir + "/FASTQ.QC.fastp.json"

	#Format conversion
	run_cmd(
		[
		"python", 
		bin + "/python/extract_bc_TELLSeq.py", 
		idxfq, 
		tmpbq
		],
		"extract_barcode"
	)
	run_cmd(
		[
		"python", 
		bin + "/python/filt_bc_TELLSeq.py", 
		"-I", tmpbq, 
		"-O1", tmpcbq, 
		"-O2", tmpebq
		],
		"extract_barcode"
	)
	
	run_cmd(
		[
		"python", 
		bin + "/python/correct_bc.TELLSeq.py", 
		"-F1", infq1, 
		"-F2", infq2, 
		"-B", idxfq, 
		"-O1", tmpfq1, 
		"-O2", tmpfq2, 
		"-C", tmpcbq
		],
		"correct_barcode"
	)
	
	#Filt & sort barcodes
	if ((sort == "Yes") & (filt == "Yes")):
		run_cmd(
			[
			"fastp", 
			"-i", tmpfq1, 
			"-I", tmpfq2, 
			"-o", tmp_filt1, 
			"-O", tmp_filt2, 
			"-w", threads, 
			"-q", "20", 
			"-u","50", 
			"-n", "10", 
			"-l", "50",
			"-h", reporthtml, 
			"-j", reportjson
			],
			"filtering_reads"
		)

		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt1 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + outfq1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmp_filt2 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + outfq2
		subprocess.call(cmd1, shell=True)
		subprocess.call(cmd2, shell=True)

		os.remove(tmpfq1)
		os.remove(tmpfq2)
		os.remove(tmp_filt1)
		os.remove(tmp_filt2)
	elif (filt == "Yes"):
		run_cmd(
			[
			"fastp", 
			"-i", tmpfq1, 
			"-I", tmpfq2, 
			"-o", outfq1, 
			"-O", outfq2, 
			"-w", threads, 
			"-q", "20", 
			"-u", "50", 
			"-n", "10", 
			"-l", "50",
			"-h", reporthtml, 
			"-j", reportjson
			],
			"filtering_reads"
		)
		os.remove(tmpfq1)
		os.remove(tmpfq2)
	elif (sort == "Yes"):
		cmd1 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmpfq1 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + outfq1
		cmd2 = "awk \'{printf(\"%s%s\",$0,(NR%4==0)?\"\\n\":\"\\t\")}\' " + tmpfq2 + " | grep \"BX:Z:\" | sort -k 2.1,2.25 -k 1,1 | tr \"\\t\" \"\\n\" > " + outfq2
		subprocess.call(cmd1, shell=True)
		subprocess.call(cmd2, shell=True)
		os.remove(tmpfq1)
		os.remove(tmpfq2)
	else:
		subprocess.call(["mv", tmpfq1, outfq1])
		subprocess.call(["mv", tmpfq2, outfq2])

	#remove tmp files
	os.remove(tmpbq)
	os.remove(tmpcbq)
	os.remove(tmpebq)


def stLFR210x(infq1, infq2, outfq1, outfq2, bin):
	######################################################################################
	# To convert the raw stLFR sequencing data to compatable 10x sequencing format       #
	######################################################################################
	if outfq1.startswith("/"):
		(outdir, filename) = os.path.split(outfq1)
	else:
		outdir = os.getcwd()

	#temp files
	tmpbq  = outdir + "/barcode_freq.txt"
	tmpmap = outdir + "/stLFRto10x.barcodes.txt"
	tmpidx = outdir + "/read-I1_si-TTCACGCG_lane-001-chunk-001.fastq"
	
	#Format conversion
	subprocess.call(["sh", bin+"/SCRIPTS/extract.bq.stLFR.sh", infq1, tmpbq])
	subprocess.call(["perl", bin+"/SCRIPTS/merge.stLFR.barcodes.pl", tmpbq, bin + "/SCRIPTS/4M-with-alts-february-2016.10x.txt", tmpmap, "8"])
	subprocess.call(["perl", bin+"/SCRIPTS/stLFR.fake.10x.pl", infq1, infq2, tmpmap, outfq1, outfq2, tmpidx])	


def TELLSeq210x(infq1, infq2, inxfq, outfq1, outfq2, white_list, bin):
	###########################################################################
	# To convert the raw TELLSeq sequencing data to 10x sequencing format     #
	###########################################################################
	if outfq1.startswith("/"):
		(outdir, filename) = os.path.split(outfq1)
	else:
		outdir = os.getcwd()

	#temp files
	tmpfq1 = outdir + "/R1_sl.fastq.gz.4tenx.fastq"
	tmpfq2 = outdir + "/R2_sl.fastq.gz.4tenx.fastq"

	#Format conversion
	subprocess.call([bin+"/TOOLS/ust10x", "-i1", inxfq, "-r1", infq1, "-r2", infq2, "-wl", white_list], cwd=outdir)
	subprocess.call(["mv", tmpfq1, outfq1])
	subprocess.call(["mv", tmpfq1, outfq1])
