import sys, getopt, json, subprocess
import numpy as np

def main(argv):
	vcffile=""
	outfile=""

	try:
		opts, args = getopt.getopt(argv,'-h:-p:-f:-o:', ["help","snvfile=","vcffile=", "outfile="])
	except getopt.GetoptError:
		print("VCF2json.py -p <input SNV vcf> -f <input sv vcf> -o <output json file>")
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == "-h":
			print("VCF2json.py -p <input SNV vcf> -f <input vcf file> -o <output json file>")
			sys.exit(2)
		elif opt in ("-p", "--snvfile"):
			snvfile=arg
		elif opt in ("-f", "--vcffile"):
			vcffile=arg
		elif opt in ("-o", "--outfile"):
			outfile=arg


	insertions={}
	deletions={}
	countsnv=0
	countlindel=0
	SNV=open(vcffile,"r")
	line=SNV.readline()
	while(line != ""):
		if(line.startswith("#")):
			pass
		else:
			LineArray = line.strip().split("\t")
			lenref = len(LineArray[3])
			lenalt = len(LineArray[4])
			if((lenref == 1) and (lenalt == 1)):
				countsnv=countsnv+1
			else:
				countlindel=countlindel+1
		line=SNV.readline()
	
	SNV.close()
	

	TotalIns=TotalDel=0
	MaxIns=MaxDel=0
	IN=open(vcffile,"r")
	record=IN.readline()
	while(record != ""):
		if(record.startswith("#")):
			pass
		else:
			LineArray = record.strip().split("\t")
			lenref = len(LineArray[3])
			lenalt = len(LineArray[4])
			lendiff = abs(lenref-lenalt) + 1
			if((lenref > lenalt) and (lenalt == 1)): ###deletion
				TotalDel=TotalDel+1
				if(lendiff in deletions.keys()):
					deletions[lendiff]=deletions[lendiff]+1
				else:
					deletions[lendiff] = 1

				if(lendiff > MaxDel):
					MaxDel = lendiff
			if((lenref < lenalt) and (lenref == 1)): ###insertion
				TotalIns=TotalIns+1
				if(lendiff in insertions.keys()):
					insertions[lendiff]=insertions[lendiff]+1
				else:
					insertions[lendiff] = 1

				if(lendiff > MaxIns):
					MaxIns = lendiff

		record=IN.readline()

	IN.close()
	
	FracINS=[0 for i in range(450)]
	FracDEL=[0 for i in range(450)]
	
	for i in range(0,450):
		length=i+50
		if(length in deletions.keys()):
			FracDEL[i]=deletions[length] / TotalDel

		if(length in insertions.keys()):
			FracINS[i]=insertions[length] / TotalIns

	VCFdata={}
	VCFdata["density_ins"]=FracINS
	VCFdata["density_del"]=FracDEL
	VCFdata["SVcount"]=TotalDel+TotalIns
	VCFdata["SNVcount"]=countsnv
	VCFdata["INDELcount"]=countlindel

	FO=open(outfile,"w")
	j=json.dumps(VCFdata, sort_keys=True, indent=2)
	FO.write(j)
	FO.close()

if __name__ == "__main__":
	main(sys.argv[1:])
