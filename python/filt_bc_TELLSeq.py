import argparse
import sys, getopt

def main():
	parser = argparse.ArgumentParser(description='Process some parameter.')
	parser.add_argument('-I',  '--infile', required=True, help='Input barcode file')
	parser.add_argument('-O1', '--outfile1', required=True, help='Output correct barcode file')
	parser.add_argument('-O2', '--outfile2', required=True, help='Output error barcode file')

	args = parser.parse_args()

	inFile  = args.infile
	outFile = args.outfile1
	errFile = args.outfile2

	count = 0
	TSBC  = open(outFile, 'w')
	ERRBC = open(errFile, 'w')
    
	IN=open(inFile,"rt")
	BarcodeDict={}
	line1=IN.readline()
	while(line1):
		line2=IN.readline()
		barcode=line2.strip()
		line3=IN.readline()
		line4=IN.readline()
		count=BarcodeDict.setdefault(barcode,0)
		BarcodeDict.update({barcode : count + 1})
		line1=IN.readline()
	
	IN.close()
	
	for key, value in BarcodeDict.items():
		if value > 1:
			TSBC.write(key + " " + str(value) +"\n")
		else:
			ErrSign=0

		for i in range(18):
			before=key[0:i]
			after=key[(i+1):18]
			if(i == 0):
				before=""
			
			if(i == 17):
				after=""
		
		seq1=before+"A"+after
		seq2=before+"C"+after
		seq3=before+"G"+after
		seq4=before+"T"+after

		if(seq1 in BarcodeDict.keys()):
			if(BarcodeDict[seq1] > 1):
				ErrSign=1

		if(seq2 in BarcodeDict.keys()):
			if(BarcodeDict[seq2] > 1):
				ErrSign=1

		if(seq3 in BarcodeDict.keys()):
			if(BarcodeDict[seq3] > 1):
				ErrSign=1

		if(seq4 in BarcodeDict.keys()):
			if(BarcodeDict[seq4] > 1):
				ErrSign=1

		if(ErrSign==1):
			ERRBC.write(key + " " + str(value) +"\n")
		else:
			TSBC.write(key + " " + str(value) +"\n")

	ERRBC.close()
	TSBC.close()



if __name__ == "__main__":
	main()

