from collections import defaultdict
import gzip
import getopt
import sys
##export genotype from vcf file
def input_vcf(in_path):
    vardict=defaultdict(list)
    f=open(in_path,"r")
    start=0
    for line in f:
        AA=line.strip('\n')
        BB=AA.split('\t')
        if BB[0]=='#CHROM':
           start=1
           continue
        if start==1 and len(BB[3])==1 and len(BB[4])==1 and BB[3]!='-' and BB[4]!='-':
           CC=BB[8].split(':')
           index=CC.index('GT')
           DD=BB[9].split(':')
           vardict[BB[0]].append([BB[1],BB[3],BB[4],DD[index]])
    f.close() 
    return vardict,list(vardict.keys())

def input_vcf_gz(in_path):
    vardict=defaultdict(list)
    f=gzip.open(in_path,"r")
    start=0
    for line in f:
        AA=line.decode().strip('\n')
        BB=AA.split('\t')
        if BB[0]=='#CHROM':
           start=1
           continue
        if start==1 and len(BB[3])==1 and len(BB[4])==1 and BB[3]!='-' and BB[4]!='-':
           CC=BB[8].split(':')
           index=CC.index('GT')
           DD=BB[9].split(':')
           vardict[BB[0]].append([BB[1],BB[3],BB[4],DD[index]])
    f.close()
    return vardict,list(vardict.keys())


def input_ref(in_path,chrlist):
    ref=defaultdict(list)
    sequence=""
    f=open(in_path,"r")
    line_index=0
    chr_index=0
    chrid=""
    for line in f:
        if line[0]!=">":
           sequence+=line.strip('\n').upper()
        else:
           if line_index!=0:
              ref[chrid]=sequence
           chrid=str(line.split()[0][1:])
           #ref[chrlist[chr_index]]=sequence
           sequence=""
        line_index+=1
           #chr_index+=1
    ref[chrid]=sequence
    f.close()
    return ref
def input_ref_gz(in_path,chrlist):
    ref=defaultdict(list)
    sequence=""
    f=gzip.open(in_path,"r")
    line_index=0
    #chr_index=0
    chrid=""
    for line in f:
        if line.decode()[0]!=">":
           sequence+=line.decode().strip('\n').upper()
        else:
           if line_index!=0:
              ref[chrid]=sequence
           chrid=str(line.split()[0][1:])
           sequence=""
           #chr_index+=1
        line_index+=1
    ref[chrid]=sequence
    f.close()
    return ref
def insertvar(ref,vcf,outprefix):
    out1=open(outprefix+'_hap1.fa','w')
    out2=open(outprefix+'_hap2.fa','w') 
    for key,value in ref.items():
        newseq1=list(value)
        newseq2=list(value)
        for pos_infor in vcf[key]:
            position=int(pos_infor[0])-1
            if pos_infor[3]=='0/1' or pos_infor[3]=='0|1' and pos_infor[3]=='1/0' or pos_infor[3]=='1|0':
               newseq1[position]=pos_infor[2]
            if pos_infor[3]=='1/1' or pos_infor[3]=='1|1':
               newseq1[position]=pos_infor[2]
               newseq2[position]=pos_infor[2]
        #newref1[key]=''.join(newseq1)
        #newref2[key]=''.join(newseq2)
        out1.write('>'+key+'\n')
        out1.write(''.join(newseq1)+'\n')
        out2.write('>'+key+'\n')
        out2.write(''.join(newseq2)+'\n')
    out1.close()
    out2.close()
    return 0
def helpinfo():
    helpinfo=\
    '''
        Generate fasta files as input of LRTK-SIM
        Version: 1.0.0
        Dependents: Python (>=3.0)
        Last Updated Date: 2017-07-22
        Contact: zhanglu2@stanford.edu

        Usage: python gen_fasta.py <options>

        Options:
                -v --vcf, the path of compressed or uncompressed vcf file
                -r --reference, the path of compressed or uncompressed ref file
                -p --prefix, prefix of new reference files
                -o --out, the path to output
                -h --help, help info
    '''
    print(helpinfo)
if __name__ == '__main__':
    if len(sys.argv) < 4:
        helpinfo()
        sys.exit(-1)
    inputvcf = None
    inputref = None
    prefix=None
    output=None

    opts, args = getopt.gnu_getopt(sys.argv[1:], 'v:r:p:o:h:', ['vcf', 'reference', 'prefix', 'out' ,'help'])
    for o, a in opts:
        if o == '-v' or o == '--vcf':
                inputvcf = a
        if o == '-r' or o == '--reference':
                inputref = a
        if o == '-p' or o == '--prefix':
                prefix = a
        if o == '-o' or o == '--out':
                output = a
        if o == '-h' or o == '--help':
                helpinfo()
                sys.exit(-1)
    if "gz" in inputvcf:
        (varlist,chrlist)=input_vcf_gz(inputvcf)
        for i in chrlist:
            print(i)
    else:
        (varlist,chrlist)=input_vcf(inputvcf)
    if "gz" in inputref:
        refseq=input_ref_gz(inputref,chrlist)
    else:
        refseq=input_ref(inputref,chrlist)
    print(refseq.keys()) 
    insertvar(refseq,varlist,output+'/'+prefix)

