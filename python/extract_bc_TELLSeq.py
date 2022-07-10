import sys
import os

L3file=sys.argv[1]
outfile=sys.argv[2]

barcode={}

IN=open(L3file,"rt")
line1=IN.readline()
while(line1):
    line2=IN.readline().strip()
    count=barcode.setdefault(line2,0)
    count = count + 1
    barcode.update({line2:count})
    line3=IN.readline()
    line4=IN.readline()
    line1=IN.readline()

IN.close()

OUT=open(outfile,"wt")
for key,value in barcode.items():
    OUT.write(key+" "+str(value)+"\n")

OUT.close()
