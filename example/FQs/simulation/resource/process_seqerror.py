#This script is used to analyze the sequencing error profile
if __name__=="__main__":
    f=open('error_profile_reads.txt',"r") # Opens file for reading
    output = open('error_profile_reads_processed.txt', 'w')
    i=0
    for line in f:
        if i==0:
            i=i+1
            output.write(line)
            continue
        clean=line.strip("\t\n")
        num=clean.split('\t')
        sequence=[]
       # print(i)
        i+=1
        for j in range(len(num)):
            if j==0 or j==1 or j==2:
               output.write(num[j])
               output.write('\t')
            elif j==6 or j==10 or j==14 or j==18:
                 sequence.append(int(num[j]))
                 sumint=sum(sequence)
                 if sumint==0:
                    for t in range(len(sequence)):
                        output.write('0')
                        output.write('\t')
                    continue
                 for t in range(len(sequence)):
                     output.write(str(sequence[t]/sumint))
                     output.write('\t')
                 sequence=[]
                 #sequence.append(int(num[j]))
            else:
           #    print(j)
               sequence.append(int(num[j]))
              # output.write(str(sequence[j]))
              # output.write('\t')
        i=i+1
        output.write('\n')
    f.close()
    output.close()
