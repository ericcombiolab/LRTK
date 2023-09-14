import sys, getopt, json, subprocess
import numpy as np

def main(argv):
	abdfile=""
	faifile=""
	vcffile=""
	outfile=""

	#print(argv)
	try:
		opts, args = getopt.getopt(argv,'-a:-f:-v:-o:', [ "abdfile=", "faifile=", "vcffile=", "outfile="])
	except getopt.GetoptError:
		print("VCF2json.py -a <input abundance file> -f <input fai file> -v <input vcf file> -o <output json file>")
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-h":
			print("VCF2json.py -a <input abundance file> -f <input fai file> -v <input vcf file> -o <output json file>")
			sys.exit(2)
		elif opt in ("-f", "--faifile"):
			faifile=arg
		elif opt in ("-v", "--vcffile"):
			vcffile=arg
		elif opt in ("-a", "--abundance"):
			abdfile=arg
		elif opt in ("-o", "--outfile"):
			outfile=arg

	ra={}
	snvs={}
	info={}

	IN=open(abdfile,"r")
	record=IN.readline()
	while(record != ""):
		if(record.startswith("#")):
			pass
		else:
			LineArray = record.strip().split("\t")
			genome=LineArray[0]
			genomeabd=LineArray[1]

			if(genome in ra.keys()):
				ra[genome]=ra[genome]+genomeabd
			else:
				ra[genome]=genomeabd

		record=IN.readline()
	IN.close()

	IN=open(faifile,"r")
	record=IN.readline()
	while(record != ""):
		if(record.startswith("#")):
			pass
		else:
			LineArray = record.strip().split("\t")
			TmpArray = LineArray[0].split("_")
			genome=TmpArray[0]+TmpArray[1]
			genomelen=LineArray[1]

			if(genome in info.keys()):
				info[genome] = info[genome] + int(genomelen)
			else:
				info[genome] = int(genomelen)

		record=IN.readline()
	IN.close()

	IN=open(vcffile,"r")
	record=IN.readline()
	while(record != ""):
		if(record.startswith("#")):
			pass
		else:
			LineArray = record.strip().split("\t")
			TmpArray = LineArray[0].split("_")
			genome=TmpArray[0]+TmpArray[1]

			if(genome in snvs.keys()):
				snvs[genome] = snvs[genome] + 1 / info[genome]
			else:
				snvs[genome] = 1 / info[genome]

		record=IN.readline()
	IN.close()

	#print(snvs)
	FracABD=[0 for i in range(450)]
	FracSNV=[0 for i in range(450)]

	#sortedABD=sorted(ra.values(), key = lambda kv:(kv[1], kv[0]))
	sortedABD=np.array(list(ra.values()))
	sortedSNV=np.array(list(snvs.values()))
	#print(sortedSNV)
	if(len(sortedABD) < 450):
		for i in range(0,len(sortedABD)):
			FracABD[i]=float(sortedABD[i])
	else:
		for i in range(0,450):
			FracABD[i]=float(sortedABD[i])

	if(len(sortedSNV) < 450):
		for i in range(0,len(sortedSNV)):
			FracSNV[i]=float(sortedSNV[i])
	else:
		for i in range(0,450):
			FracSNV[i]=float(sortedSNV[i])

	FracABD.sort(reverse=True)
	FracSNV.sort(reverse=True)

	VCFdata={}
	VCFdata["abd"]=FracABD
	VCFdata["snv"]=FracSNV

	FO=open(outfile,"w")
	j=json.dumps(VCFdata, sort_keys=True, indent=2)
	FO.write(j)
	FO.close()

if __name__ == "__main__":
	main(sys.argv[1:])
