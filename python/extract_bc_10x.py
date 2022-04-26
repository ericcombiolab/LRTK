
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
            tmp = count // 4
            if(tmp % 2 == 0):
                tmp = count + 1
                if((tmp % 2) == 0):
                    bc = line[:16] + '\n'
                    tenX_bc.write(bc)
                    bases = line[23:]
                    tenX.write(bases)
                else:
                    tenX_bc.write(line)
                    tenX.write(line)
            else:
                tenX.write(line)
            count += 1

    tenX_bc.close()
    tenX.close()


if __name__ == "__main__":
#    print(sys.argv[1:])
    main(sys.argv[1:])

