import os, sys, getopt, json, subprocess
import numpy as np

def main(argv):
	fragment_file=""
	stat_file=""
	outfile=""
	
	try:
		opts, args = getopt.getopt(argv,'-h:-f:-s:-o:', ["help","fragment=", "stat=", "outfile="])
	except getopt.GetoptError:
		print("Fragment2json.py -f <input fragment file> -s <input statistic file> -o <output json file>")
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == "-h":
			print("Fragment2json.py -f <input fragment file> -s <input statistic file> -o <output json file>")
			sys.exit()
		elif opt in ("-f", "--fragment"):
			fragment_file=arg
		elif opt in ("-s", "--stat"):
			stat_file=arg
		elif opt in ("-o", "--outfile"):
			outfile=arg

	if outfile.startswith("/"):	
		(outdir, filename) = os.path.split(outfile)
	else:
		outdir = os.getcwd()
		(outdir, filename) = os.path.split(outdir + "/" + outfile)

	BamData={}
	ST=open(stat_file,"r")
	line=ST.readline()
	while(line != ""):
		line=ST.readline()
		if(line.find("alignments") > 0):
			next
		if(line.find(":") > 0):
			LineArray = line.strip().split(":")
			if(LineArray[0] == "Total number of read pairs"):
				BamData["TotalReads"]=LineArray[1].strip()
			elif(LineArray[0] == "Duplication ratio"):
				TmpArray = LineArray[1].strip().split("=")
				BamData["DuplicationRate"]=TmpArray[1].strip()
			elif(LineArray[0] == "Total number of barcodes"):
				BamData["TotalBarcodes"]=LineArray[1].strip()
			elif(LineArray[0] == "Total number of fragments"):
				BamData["TotalFragments"]=LineArray[1].strip()
			elif(LineArray[0] == "Fragment coverage (C_F)"):
				BamData["C_F"]=LineArray[1].strip()
			elif(LineArray[0] == "Read coverage per fragment (C_R)"):
				BamData["C_R"]=LineArray[1].strip()
			elif(LineArray[0] == "Read coverage (C)"):
				BamData["C"]=LineArray[1].strip()
			elif(LineArray[0] == "Mean fragment length (U_FL)"):
				BamData["U_FL"]=LineArray[1].strip()
			elif(LineArray[0] == "Fragment length max, N25, N50 (WU_FL), N75, min"):
				TmpArray = LineArray[1].strip().split(" ")
				BamData["WU_FL"]=TmpArray[1].strip().replace(",","")
			elif(LineArray[0] == "Fragments per barcode (N_F/P)"):
				BamData["NF_P"]=LineArray[1].strip()
			elif(LineArray[0] == "Median insert size"):
				BamData["InsertSize"]=LineArray[1].strip()
	
	ST.close()

	###Processing Fragment file
	nfp_barcode_last   = "";
	nfp_fragment_count = 0;
	nfp_fragment_total = 0;
	nfp_barcode_total  = 0;
	nfp={};
	###end

	###Unweighted and Weighted Mean Fragment Length
	ufl_FragmentLength_total= 0
	ufl_FragmentCount_total= 0
	ufl_FragmentLength = {}
	ufl_FragmentLength_Max = 0
	###end

	###Read Coverage of fragment
	cr_FragmentLength_total = 0;
	cr_ReadCount_total = 0;
	cr_ReadCoverage_PerFragment = 0;
	cr_ReadCoverage={};
	cr_FragmentCount_total = 0;
	###end

	###Fragment Coverage
	fr_FragmentLength_total = 0;
	fr_GenomicSpace={}
	fr_FragmentCoverage=0;
	###end
	
	sorted_fragment_file = fragment_file + ".sorted"
	cmd = "sort -k2 -k3 -k4 " + fragment_file + " > " + sorted_fragment_file
	subprocess.call(cmd, shell=True)
	FR=open(sorted_fragment_file,"r")
	record=FR.readline()
	while(record != ""):
		if(record.startswith("ID")):
			pass
		else:
			LineArray = record.strip().split()

			###NFP start
			nfp_fragment_total = nfp_fragment_total +1
			nfp_barcode_current = LineArray[1]
			if(nfp_barcode_current == nfp_barcode_last):
				nfp_fragment_count = nfp_fragment_count + 1
			else:
				if(nfp_barcode_last == ""):
					nfp_barcode_total = nfp_barcode_total + 1
					nfp_barcode_last = nfp_barcode_current
				else:
					if(nfp_fragment_count in nfp.keys()):
						nfp[nfp_fragment_count] = nfp[nfp_fragment_count] + 1
					else:
						nfp[nfp_fragment_count] = 1

					nfp_barcode_total = nfp_barcode_total + 1
					nfp_fragment_count = 0
					nfp_barcode_last = nfp_barcode_current
			###NFP end
			###uFL start
			ufl_FragmentCount_total = ufl_FragmentCount_total + 1
			ufl_FragmentLength_short = int(int(LineArray[5]) / 1000)
			if(ufl_FragmentLength_short in ufl_FragmentLength.keys()):
				ufl_FragmentLength[ufl_FragmentLength_short] = ufl_FragmentLength[ufl_FragmentLength_short] + 1
			else:
				ufl_FragmentLength[ufl_FragmentLength_short] = 1

			if(ufl_FragmentLength_short > ufl_FragmentLength_Max):
				ufl_FragmentLength_Max = ufl_FragmentLength_short

	        	###uFL end
			###cr start
			cr_FragmentCount_total = cr_FragmentCount_total + 1
			cr_coverage_rate = round(int(LineArray[6]) * 150 / int(LineArray[5]),1)
			if(cr_coverage_rate in cr_ReadCoverage.keys()):
				cr_ReadCoverage[cr_coverage_rate] = cr_ReadCoverage[cr_coverage_rate] + 1
			else:
				cr_ReadCoverage[cr_coverage_rate] = 1;
			###cr end

			###fr start
			###fr end
		#print(record.strip())
		record=FR.readline()

	FR.close()

	NFPfile = outdir + "/Density.NFP.xls"
	OUTNFP  = open(NFPfile,"w")
	for i in range(20):
		number=i+1
		if(number in nfp.keys()):
			nfp_value = nfp[i] / nfp_barcode_total
		else:
			nfp_value = 0

		OUTNFP.write("\t".join([str(number),str(nfp_value),"\n"]))
	
	OUTNFP.close()

	FLfile = outdir + "/Density.FL.xls"
	OUTFL=open(FLfile, "w")
	for i in range(400):
		number =i+1
		if(number in ufl_FragmentLength.keys()):
			ufl_value = ufl_FragmentLength[number] / ufl_FragmentCount_total
		else:
			ufl_value = 0
		
		OUTFL.write("\t".join([str(number),str(ufl_value),"\n"]))

	OUTFL.close()

	CRfile = outdir + "/Density.CR.xls"
	OUTCR=open(CRfile,"w")
	for i in range(100):
		number=(i+1)/10
		if(number in cr_ReadCoverage.keys()):
			cr_value = cr_ReadCoverage[number] / cr_FragmentCount_total
		else:
			cr_value = 0
		
		OUTCR.write("\t".join([str(number),str(cr_value),"\n"]))

	OUTCR.close()

	FO=open(outfile,"w")
	j=json.dumps(BamData, sort_keys=True, indent=2)
	FO.write(j)
	FO.close()


if __name__ == "__main__":
	main(sys.argv[1:])
