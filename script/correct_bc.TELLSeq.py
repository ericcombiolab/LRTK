#!/usr/bin/python
import argparse
import os
import gzip
import sys

def main():
	parser = argparse.ArgumentParser(description='Process some parameter.')
	parser.add_argument('-F1', '--fastq1',   required=True, help='Input FASTQ file')
	parser.add_argument('-F2', '--fastq2',   required=True, help='Input FASTQ file')
	parser.add_argument('-B',  '--barcode',  required=True, help='Input BARCODE file')
	parser.add_argument('-O1', '--outfile1', required=True, help='Output correct FASTQ file')
	parser.add_argument('-O2', '--outfile2', required=True, help='Output correct FASTQ file')
	parser.add_argument('-C',  '--correct',  required=True, help='Correct barcodes')

	args = parser.parse_args()
	barcodefile=args.barcode
	FQ1file=args.fastq1
	FQ2file=args.fastq2
	outfq1=args.outfile1
	outfq2=args.outfile2
	dbfile=args.correct

	IN=open(barcodefile,"rt")
	FI1=open(FQ1file,"rt")
	FI2=open(FQ2file,"rt")
	FO1=open(outfq1,"wt")
	FO2=open(outfq2,"wt")
	DB=open(dbfile,"rt")

	barcodes={}
	line1=DB.readline()
	while(line1):
		array=line1.strip().split(" ")
		barcodes.update({array[0]:1})
		line1=DB.readline()
	
	DB.close()

	if 'AAAAAAAAAAAAAAAAAA' in barcodes.keys():
		log=barcodes.pop('AAAAAAAAAAAAAAAAAA')

	if 'CCCCCCCCCCCCCCCCCC' in barcodes.keys():
		log=barcodes.pop('CCCCCCCCCCCCCCCCCC')

	if 'GGGGGGGGGGGGGGGGGG' in barcodes.keys():
		log=barcodes.pop('GGGGGGGGGGGGGGGGGG')

	if 'TTTTTTTTTTTTTTTTTT' in barcodes.keys():
		log=barcodes.pop('TTTTTTTTTTTTTTTTTT')


	print("Processing FQs...")
	bline1=IN.readline()
	while(bline1):
		readid=bline1.strip().split(" ")
		bline2=IN.readline().strip()
		bline3=IN.readline()
		bline4=IN.readline()
		fq1line1=FI1.readline()
		fq1readid=fq1line1.strip().split(" ")
		newfq1line1=fq1readid[0]+" BX:Z:"+bline2+"-1\n"
		fq1line2=FI1.readline()
		fq1line3=FI1.readline()
		fq1line4=FI1.readline()
		fq2line1=FI2.readline()
		fq2readid=fq2line1.strip().split(" ")
		newfq2line1=fq2readid[0]+" BX:Z:"+bline2+"-1\n"
		fq2line2=FI2.readline()
		fq2line3=FI2.readline()
		fq2line4=FI2.readline()

		if bline2 in barcodes.keys():
			FO1.write(newfq1line1)
			FO1.write(fq1line2)
			FO1.write(fq1line3)
			FO1.write(fq1line4)
			FO2.write(newfq2line1)
			FO2.write(fq2line2)
			FO2.write(fq2line3)
			FO2.write(fq2line4)
		else:
			FO1.write(fq1readid[0])
			FO1.write("\n")
			FO1.write(fq1line2)
			FO1.write(fq1line3)
			FO1.write(fq1line4)
			FO2.write(fq2readid[0])
			FO2.write("\n")
			FO2.write(fq2line2)
			FO2.write(fq2line3)
			FO2.write(fq2line4)

		bline1=IN.readline()

	IN.close()
	FI1.close()
	FI2.close()
	FO1.close()
	FO2.close()

if __name__ == "__main__":
	main()

