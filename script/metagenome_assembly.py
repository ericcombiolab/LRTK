import os
import sys
import subprocess
import shutil

python_path = os.path.dirname(os.path.abspath( __file__ )) + "/"
sys.path.append(python_path)

pangaea_path = python_path + "Pangaea/"
from utility import *

def MAGconstruction(fq1, fq2, metaspades, athena_local, athena_hybrid, lt, cluster, outfile, threads):
	###########################################################################
	# To align the clean ULRF sequencing reads to reference genome            #
	###########################################################################
	if outfile.startswith("/"):
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir,filename)=os.path.split(outdir + "/" + outfile)
	
	#temp paths and files
	assembled_genome  = outdir + "/final.asm.fa"
	contig = outdir + "/uncircularized.fa"
	circular_dir = outdir + "/4.assembly/quickmerge/circular/3.circularization/3.circular_sequences/"

	res_binning  = outdir + "/binning"
	os.makedirs(res_binning)

	tmp_ref      = res_binning + "/reference"
	tmp_sam      = res_binning + "/tmp.sam"
	tmp_bam      = res_binning + "/tmp.bam"
	tmp_depth    = res_binning + "/final.depth.txt"
	
	#assembly reads
	cmd = "python " + pangaea_path + "/pangaea.py -1 " + fq1 + " -2 " + fq2 + " -sp " + metaspades + " -lc " + athena_local + " -at " + athena_hybrid + " -lt 10,30 -c 30 -o " + outdir
	#print(cmd)
	run_cmd(cmd, "pangaea")
	
	#preprocess
	hash_dict = {}
	files = [f for f in os.listdir(indir) if f.endswith('.fa')]
	for contigfile in files:
		name = os.path.basename(contigfile)
		tmp = name.split('_')
		contig_id = ">{}_{}".format(tmp[0], tmp[1])
		hash_dict[contig_id] = 1

	db = {}
	with open(fafile, 'r') as infile, open(outfile, 'w') as outfile:
		contig_id = infile.readline().strip()
		line = infile.readline().strip()
		contig_seq = "NA"

		while line:
			line = line.strip()

			if line.startswith('>'):
				outfile.write(contig_id + "\n")
				outfile.write(contig_seq + "\n")
				contig_seq = "NA"
				contig_id = line
			else:
				if contig_seq == "NA":
					contig_seq = line
				else:
					contig_seq = contig_seq + "\n" + line

		line = infile.readline()

	outfile.write(contig_seq + "\n")
	infile.close()
	outfile.close()

	#binning
	cmd = "bowtie2-build -f " + contig + " " + tmp_ref + " --threads " + str(threads)
	#print(cmd)
	subprocess.call(cmd, shell=True)

	cmd = "bowtie2 -1 " + fq1 + " -2 " + fq2 + " -p " + str(threads) + " -x " + tmp_ref + " -S " + tmp_sam
	#print(cmd)
	subprocess.call(cmd, shell=True)

	cmd = "samtools sort -@ " + str(threads) + " -l 9 " + " -O BAM " + tmp_sam + " -o " + tmp_bam
	#print(cmd)
	subprocess.call(cmd, shell=True)

	cmd = "jgi_summarize_bam_contig_depths --outputDepth " + tmp_depth + " " + tmp_bam
	#print(cmd)
	subprocess.call(cmd, shell=True)

	cmd = "metabat2 -m 1500 -t " + str(threads) + " -i " + tmp_ref + " -a " + tmp_depth + " -o " + res_binning + "/all -v"
	print(cmd)
	run_cmd(cmd, "metabat2")

	#remove tmp files
	os.remove(tmp_ref)
	os.remove(tmp_sam)
	os.remove(tmp_bam)
	os.remove(tmp_depth)
