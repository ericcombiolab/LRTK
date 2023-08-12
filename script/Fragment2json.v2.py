import sys, getopt, json, subprocess
import numpy as np

def main(argv):
	bam_file=""
	fragment_file=""
	insert_file=""
	stat_file=""
	samtools=""
	outfile=""
	
	try:
		opts, args = getopt.getopt(argv,'-h:-b:-f:-i:-s:-t:-o:', ["help","bam_file=","fragment=","insert_file=", "stat=", "samtools=","outfile="])
	except getopt.GetoptError:
		print("Fragment2json.py -b <input alignment file> -f <input fragment file> -i <insert file> -s <input statistic file> -t <path to samtools> -o <output json file>")
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == "-h":
			print("Fragment2json.py -b <input alignment file> -f <input fragment file> -i <insert file> -s <input statistic file> -t <path to samtools> -o <output json file>")
			sys.exit()
		elif opt in ("-b", "--bam_file"):
			bam_file=arg
		elif opt in ("-f", "--fragment"):
			fragment_file=arg
		elif opt in ("-i", "--insert_file"):
			insert_file=arg
		elif opt in ("-s", "--stat"):
			stat_file=arg
		elif opt in ("-t", "--samtools"):
			samtools=arg
		elif opt in ("-o", "--outfile"):
			outfile=arg

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
#		print(record.strip())
		record=FR.readline()

	FR.close()

	nfp_array = [0 for i in range(50)]
	for i in range(50):
		if(i in nfp.keys()):
			nfp_array[i] = nfp[i] / nfp_barcode_total
	
	BamData["density_nfp"]=nfp_array

	fl_array = [0 for i in range(100)]
	for i in range(100):
		if(i in ufl_FragmentLength.keys()):
			fl_array[i] = ufl_FragmentLength[i] / ufl_FragmentCount_total
        
	BamData["density_fl"]=fl_array

	cr_array = [0 for i in range(100)]
	for i in range(100):
		val=i/10
		if(val in cr_ReadCoverage.keys()):
			cr_array[i]= cr_ReadCoverage[val] / cr_FragmentCount_total

	BamData["density_cr"]=cr_array

	print("processing depth ...")
	depthfile = outfile + ".depth.xls"
	cmd = samtools + " depth -q 20 -Q 20 " + bam_file + " > " + depthfile
	subprocess.call(cmd, shell=True)
	COVdata={}
	TotalPOS=0
	DEP=open(depthfile,"r")
	record=DEP.readline()
	while(record != ""):
		LineArray = record.strip().split("\t")
		TotalPOS=TotalPOS+1
		if(LineArray[2] in COVdata.keys()):
			value = LineArray[2].strip()
			COVdata[value] = COVdata[value] + 1
		else:
			value = LineArray[2].strip()
			COVdata[value] = 1

		record=DEP.readline()

	DEP.close()

	#print(COVdata.keys())
	COVarray = [0 for i in range(200)]
	for i in range(200):
		number = str(i +1)
		if(number in COVdata.keys()):
			#print(COVdata[number])
			COVarray[i] = COVdata[number] / TotalPOS
	
	BamData["depth"] = COVarray

	print("processing insert size ...")
	INS=open(insert_file,"r")
	INSdata={}
	Total_ins_record=0
	record=INS.readline()
	while(record != ""):
		if(record.startswith("Insert_Size")):
			pass
		else:
			Total_ins_record=Total_ins_record+1
			size=record.strip()
			if(size in INSdata.keys()):
				INSdata[size]=INSdata[size]+1
			else:
				INSdata[size]=1

		record=INS.readline()

	INS.close()
	
	INSarray = [0 for i in range(700)]
	for i in range(700):
		number = str(i +1)
		if(number in INSdata.keys()):
			#print(INSdata[number])
			INSarray[i] = INSdata[number] / Total_ins_record

	BamData["insert_size"] = INSarray

	FO=open(outfile,"w")
	j=json.dumps(BamData, sort_keys=True, indent=2)
	FO.write(j)
	FO.close()


if __name__ == "__main__":
	main(sys.argv[1:])
