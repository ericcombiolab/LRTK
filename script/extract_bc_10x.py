
import sys, getopt

def main(argv):
    inFile =''
    outFile =''
    outBarcode = ''
    try:
        opts, args = getopt.getopt(argv, "", ["help", "infile=", "outfile=", "outbarcode="])
    except getopt.GetoptError:
        print("extract_bc.py -i <input fastq file> -o <output fastq file> -b <output barcode file>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("extract_bc.py -i <input fastq file> -o <output fastq file> -b <output barcode file>")
            sys.exit()
        elif opt in ("-i", "--infile"):
            inFile = arg
        elif opt in ("-o", "--outfile"):
            outFile = arg
        elif opt in ("-b", "--outbarcode"):
            outBarcode = arg

    count = 0
    tenX_bc = open(outBarcode, 'w')
    tenX = open(outFile, 'w')
    
    with open(inFile,'rt') as reads:
        for line in reads:
            if(count == 0 or count == 2):
                tenX_bc.write(line)
                tenX.write(line)
                count=count+1
            elif(count == 1):
                bc = line[:16] + '\n'
                tenX_bc.write(bc)
                bases = line[23:]
                tenX.write(bases)
                count=count+1
            elif(count == 3):
                bc = line[:16] + '\n'
                tenX_bc.write(bc)
                bases = line[23:]
                tenX.write(bases)
                count=0

    tenX_bc.close()
    tenX.close()


if __name__ == "__main__":
#    print(sys.argv[1:])
    main(sys.argv[1:])

