from argparse import ArgumentParser
import time
import os

parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n liuyichen@std.uestc.edu.cn\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--ref_fa','-r',help="Reference fasta file for reformat",required=True)
parser.add_argument('--in_vcf','-i',help="Original vcf file",required=True)
parser.add_argument('--out_vcf','-o',help="Output reformated vcf file",required=True)
parser.add_argument('--add_header','-head',help="Add header to vcf (38 or 19)",default=False)
parser.add_argument('--add_chr','-ac',help="Add 'chr' to CHROM field (1->chr1)",action="store_true")
parser.add_argument('--gz_tbi','-gt',help="Output gz and tbi file",action="store_true")
parser.add_argument('--base_norm','-bn',help="If set, change base from 0 to 1 (for all types of variants) and add 1 base at the beginning of INDEL and SV.",action="store_true")
args = parser.parse_args()

script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 


def GetGenoSeq (fafile):
    genome ={}
    chro = []
    chrnum = ""
    with open(fafile,"r") as f:
        for line in f:
            if line.startswith('>'):
                if chro:
                    genome[chrnum] = ''.join(chro)
                    chro = []
                chrnum = line[1:].split(' ')[0].strip('\n')
            else:
                chro.append(line.strip('\n'))
        genome[chrnum] = ''.join(chro)
        chro = []
    return genome
    
def modify(vcf_file,genome,vcfwrite,add_chr,add_header,pos_base):
    chrom_list = ['chr1','chr2','chr3','chr4','chr5',
                  'chr6','chr7','chr8','chr9','chr10',
                  'chr11','chr12','chr13','chr14','chr15',
                  'chr16','chr17','chr18','chr19','chr20',
                  'chr21','chr22','chrX']
    with open (vcf_file, "r") as f:
        with open (vcfwrite,"w") as fw:
            if add_header:
                with open(code_path+"header/header"+add_header,"r") as fh:
                    for line in fh:
                        fw.write(line)
            for line in f:
                if line[0] == '#':
                    fw.write(line)
                else:
                    line = line.split('\t')
                    POS = int(line[1])
                    if add_chr:
                        line[0] = "chr" + line[0]
                    CHROM = line[0]
                    INFO = line[7]
                    line[3] = line[3].upper()
                    line[4] = line[4].upper()
                    if CHROM in chrom_list:
                        if INFO == 'SVTYPE=SNP':
                            if pos_base:
                                line[1] = str(POS+1)
                            fw.write('\t'.join(line))
                        elif INFO == 'SVTYPE=DEL':
                            if pos_base:
                                line[3] = genome[CHROM][POS-1].upper() + line[3]
                                line[4] = genome[CHROM][POS-1].upper()
                                line[1] = str(POS)
                            fw.write('\t'.join(line))
                        else:#INFO == 'VTYPE=INS'
                            if pos_base:
                                line[3] = genome[CHROM][POS-1].upper()
                                line[4] = genome[CHROM][POS-1].upper() + line[4]
                                line[1] = str(POS)
                            fw.write('\t'.join(line))

def GzTbi(vcfwrite):
    vcfsort = os.popen('vcf-sort '+vcfwrite+' > '+vcfwrite[:-4]+'_sorted.vcf')
    print(vcfsort.read())
    bgzip = os.popen('bgzip -c '+vcfwrite[:-4]+'_sorted.vcf'+' > '+vcfwrite[:-4]+'_sorted.vcf.gz')
    print(bgzip.read())
    tabix = os.popen('tabix -p vcf '+vcfwrite[:-4]+'_sorted.vcf.gz')
    print(tabix.read())

if __name__ == "__main__":
    ref_genome = args.ref_fa
    vcf_file = args.in_vcf
    vcfwrite = args.out_vcf
    add_chr = args.add_chr
    add_header = args.add_header
    gz_tbi = args.gz_tbi
    pos_base = args.base_norm
    print("Vcf reformat start")
    t = time.time()
    genome = GetGenoSeq(ref_genome)
    modify(vcf_file,genome,vcfwrite,add_chr,add_header,pos_base)
    if gz_tbi:
        GzTbi(vcfwrite)
    print("Vcf reformat finished")
    print("Time used:",time.time()-t)
