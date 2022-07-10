import gzip
import sys

def main(argv):
    inFile = argv[1]
    outFile = argv[2] + "_read2_only.fq"
    outbc1 = argv[2] + "_bc1.fq"
    outbc2 = argv[2] + "_bc2.fq"
    outbc3 = argv[2] + "_bc3.fq"

    bc1 = open(outbc1, 'w')
    bc2 = open(outbc2, 'w')
    bc3 = open(outbc3, 'w')

    lfr_read2 = open(outFile, 'w')
    count = 0
    with open(inFile,'rt') as read2:

        for line in read2:
            count +=1 
            if(count % 2 == 0):
                b1 = line[100:110] + '\n'
                b2 = line[116:126] + '\n'
                b3 = line[132:] + '\n'
                bc1.write(b1)
                bc2.write(b2)
                bc3.write(b3)
                bases = line[:100] + '\n' 
                lfr_read2.write(bases)
            else:
                bc1.write(line)
                bc2.write(line)
                bc3.write(line)
                lfr_read2.write(line)
    bc1.close()
    bc2.close()
    bc3.close()
    read2.close()	

if __name__ == "__main__":
    main(sys.argv)
