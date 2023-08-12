import sys, getopt, json, subprocess
import numpy as np

def main(argv):
	abdfile=""
	vcffile=""
	outfile=""
    
    print(argv)
	try:
		opts, args = getopt.getopt(argv,'-a:-f:-o:', [ "abdfile=", "vcffile=", "outfile="])
	except getopt.GetoptError:
		print("VCF2json.py -a <input abundance file> -f <input vcf file> -o <output json file>")
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == "-h":
			print("VCF2json.py -a <input abundance file> -f <input vcf file> -o <output json file>")
			sys.exit(2)
		elif opt in ("-f", "--vcffile"):
			vcffile=arg
		elif opt in ("-a", "--abundance"):
			abdfile=arg
		elif opt in ("-o", "--outfile"):
			outfile=arg

	ra={}
	snvs={}
        
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

		IN.close()

	FracABD=[0 for i in range(450)]
	FracDEL=[0 for i in range(450)]
	sortedABD=sorted(ra.items(), key = lambda kv:(kv[1], kv[0]))
	for i in range(0,450):
		FracABD[i]=sortedABD[i]

	VCFdata={}
	VCFdata["abd"]=FraABD

	FO=open(outfile,"w")
	j=json.dumps(VCFdata, sort_keys=True, indent=2)
	FO.write(j)
	FO.close()

if __name__ == "__main__":
	main(sys.argv[1:])
