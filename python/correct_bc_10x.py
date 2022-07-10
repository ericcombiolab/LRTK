import sys, getopt

def main(argv):
        WHITELIST=''
        SAMFILE=''
        READ1=''
        READ2=''
        OUTFQ1=''
        OUTFQ2=''
        NM=''

        try:
                opts, args = getopt.getopt(argv,'-h-w:-s:-1:-2:-3:-4:-n:' ,["help", "whitelist=", "samfile=", "read1=", "read2=", "outfq1=", "outfq2=", "nm="])
        except getopt.GetoptError:
                print("correct_bc_10x.py -1 <input fastq file 1> -2 <input fastq file 2> -s <input sam file>  -w <white_list_10x_barcode> -3 <output fastq file 1> -4 <output fastq file 2> -n <number of mismatches>")
                sys.exit(2)

        for opt, arg in opts:
                if opt == "-h":
                        print("correct_bc_10x.py -1 <input fastq file 1> -2 <input fastq file 2> -s <input sam file> -w <white_list_10x_barcode> -3 <output fastq file 1> -4 <output fastq file 2> -n <number of mismatches>")
                        sys.exit()
                elif opt in ("-w", "--whitelist"):
                        WHITELIST=arg
                elif opt in ("-s", "--samfile"):
                        SAMFILE=arg
                elif opt in ("-1", "--read1"):
                        READ1=arg
                elif opt in ("-2", "--read2"):
                        READ2=arg
                elif opt in ("-3", "--outfq1"):
                        OUTFQ1=arg
                elif opt in ("-4", "--outfq2"):
                        OUTFQ2=arg
                elif opt in ("-n", "--nm"):
                        NM=arg

        ref={}
        ref_alter={}
        RI=open(WHITELIST,"r")
        line_id=RI.readline()
        ref_id=0
        while(line_id):
                if line_id.startswith(">"):
                        ref_id=line_id.strip().replace(">","")
                line_barcode=RI.readline()
                line_barcode=line_barcode.strip()
                ref.update({line_barcode:ref_id})
                ref_alter.update({ref_id:line_barcode})
                line_id=RI.readline()
        
        RI.close()

        SI=open(SAMFILE,"r")    
        FI1=open(READ1,"r")
        FI2=open(READ2,"r")
        FO1=open(OUTFQ1,"w")
        FO2=open(OUTFQ2,"w")
        count_match=count_filter=0
        type_barcode={}

        line=SI.readline()
        while(line):
                if not line.startswith("@"): # ignoring the header in samfile
                        fq1_line1=FI1.readline()
                        fq1_line2=FI1.readline()
                        fq1_line3=FI1.readline()
                        fq1_line4=FI1.readline()
                        fq2_line1=FI2.readline()
                        fq2_line2=FI2.readline()
                        fq2_line3=FI2.readline()
                        fq2_line4=FI2.readline()
                        
                        array = line.strip().split("\t") #
                        read_id=array[0]
                        barcode=array[9]
                        barcode_com_rev_array=["#" for i in range(16)]
                        for i in range(16):
                                if barcode[i]=="A":
                                        barcode_com_rev_array[15-i]="T"
                                elif barcode[i]=="C":
                                        barcode_com_rev_array[15-i]="G"
                                elif barcode[i]=="G":
                                        barcode_com_rev_array[15-i]="C"
                                elif barcode[i]=="T":
                                        barcode_com_rev_array[15-i]="A"

                        barcode_com_rev="".join(barcode_com_rev_array)
                        flag=bin(int(array[1]))
						
                        if len(array) < 12: ##status: completely unmapped
                                read_id=array[0]
                        elif ((barcode in ref.keys()) or (barcode_com_rev in ref.keys())): ##barcode in whitelist:soft clipping reads
                                if barcode in type_barcode.keys():
                                        type_barcode.update({barcode:(type_barcode[barcode]+1)})
                                else:
                                        type_barcode.setdefault(barcode, 1)

                                barcode="".join(["BX:Z:",barcode,"-1"])
                        elif ((len(flag) < 5) or (int(flag[-3]) == 0)):  ##Correctly mappped  and NM <= 2
                                NMfield=array[12].split(":")
                                if (int(NMfield[2]) <= int(NM)):
                                        barcode="".join(["BX:Z:",ref_alter[array[2]],"-1"])
                                if ref_alter[array[2]] in type_barcode.keys():
                                        type_barcode.update({ref_alter[array[2]]:(type_barcode[ref_alter[array[2]]]+1)})
                                else:
                                        type_barcode.setdefault(ref_alter[array[2]], 1)
                                
                        if barcode.startswith("BX"):
                                count_match=count_match+1
                                read_id="@"+array[0]+" "+barcode+"\n"
                        else:
                                count_filter=count_filter+1
                                read_id="@"+array[0]+"\n"

                        FO1.write(read_id)
                        FO1.write(fq1_line2)
                        FO1.write(fq1_line3)
                        FO1.write(fq1_line4)
                        
                        FO2.write(read_id)
                        FO2.write(fq2_line2)
                        FO2.write(fq2_line3)
                        FO2.write(fq2_line4)

                line=SI.readline()
        SI.close()
        FI1.close()
        FI2.close()
        FO1.close()
        FO2.close()

#        print("#SUMMARY:MatchedBarcodes: %d" % count_match)
#        print("#SUMMARY:FilteredBarcodes: %d" % count_filter)
#        print("#SUMMARY:BARCODES\tCount")
#        for key,value in type_barcode.items():
#                print("%s\t%d" % (key, value))

if __name__ == "__main__":
        main(sys.argv[1:])
