#include <htslib/sam.h>
#include <iostream> 
#include <stdio.h> 
#include <string.h> 
#include <string> 
#include <map> 
#include <fstream>
#include <cmath>
using namespace std;

int main(int arg, char *argvs[])
{
    //argvs[1] barcode whitelist bc.fa
    //argvs[2] aligned barcode1 bc2.aln.sam
    //argvs[3] aligned barcode2 bc2.aln.sam
    //argvs[4] aligned barcode3 bc3.aln.sam
    //argvs[5] 1st reads without barcode read1.fq
    //argvs[6] 2nd reads withtou barcode read2.fq
    //argvs[7] number of mismatch allowed NM
    //argvs[8] output1 file name for the read 1
    //argvs[9] output1 file name for the read 2
    
    printf("Get barcode whitelistes\n");
    map <string, string> ref_bc;  
    fstream barcode;
    barcode.open(argvs[1], ios::in);
    string temp;
    while(getline(barcode, temp))
    {
        string refName = temp.substr(1, temp.length());
        getline(barcode, temp);
        ref_bc.insert({refName, temp});
    }

    samFile *fp1 = hts_open(argvs[2], "r");
    samFile *fp2 = hts_open(argvs[3], "r");
    samFile *fp3 = hts_open(argvs[4], "r");
    //int count = 0; 
    bam_hdr_t *bamHdr1 = sam_hdr_read(fp1);
    bam1_t *aln1 = bam_init1();
    bam_hdr_t *bamHdr2 = sam_hdr_read(fp2);
    bam1_t *aln2 = bam_init1();
    bam_hdr_t *bamHdr3 = sam_hdr_read(fp3);
    bam1_t *aln3 = bam_init1();

    fstream output1;
    output1.open(argvs[8], ios::out);
    fstream output2;
    output2.open(argvs[9], ios::out);

    fstream reads1;
    reads1.open(argvs[5], ios::in);
    fstream reads2;
    reads2.open(argvs[6], ios::in);
    string temp1;

    while(sam_read1(fp1, bamHdr1, aln1) == 0 && sam_read1(fp2, bamHdr2, aln2) == 0 && sam_read1(fp3, bamHdr3, aln3) == 0 && getline(reads1, temp) && getline(reads2, temp1))
    {
        string seq;
        getline(reads1, seq);
        string sym;
        getline(reads1, sym);
        string qual;
        getline(reads1, qual);

        string seq1;
        getline(reads2, seq1);
        string sym1;
        getline(reads2, sym1);
        string qual1;
        getline(reads2, qual1);
    
        char *flag1 = bamHdr1->target_name[aln1->core.flag];
        char *flag2 = bamHdr2->target_name[aln2->core.flag];
        char *flag3 = bamHdr3->target_name[aln3->core.flag];

        if(strcmp(flag1, "5") == 0 || flag1 == NULL || strcmp(flag2, "5") == 0 || flag2 == NULL || strcmp(flag3, "5") == 0 || flag3 == NULL)     
        {
            string line1 = temp;
            output1 << line1 << endl;
            output1 << seq << endl;
            output1 << sym << endl;
            output1 << qual << endl;        


            string line11 = temp1;
            output2 << line11 << endl;
            output2 << seq1 << endl;
            output2 << sym1 << endl;
            output2 << qual1 << endl;
            continue;
        }

        char *rname1 = bamHdr1->target_name[aln1->core.tid];
        char *rname2 = bamHdr2->target_name[aln2->core.tid];
        char *rname3 = bamHdr3->target_name[aln3->core.tid];

        uint8_t *NM1 = bam_aux_get(aln1, "NM");
        int64_t nm1 = bam_aux2i(NM1);
        uint8_t *NM2 = bam_aux_get(aln2, "NM");
        int64_t nm2 = bam_aux2i(NM2);
        uint8_t *NM3 = bam_aux_get(aln3, "NM");
        int64_t nm3 = bam_aux2i(NM3);
        
        int NM = atoi(argvs[7]);
        if(nm1 + nm2 + nm3 <= NM)
        {
            auto bc1 = ref_bc.find(rname1);
            auto bc2 = ref_bc.find(rname2);
            auto bc3 = ref_bc.find(rname3);

            if(bc1 != ref_bc.end() && bc2 != ref_bc.end() && bc3 != ref_bc.end())
            {
                string seq_bc = " BX:Z:" + bc1->second + bc2->second + bc3->second;
                string line1 = temp + seq_bc + "-1";
                output1 << line1 << endl;
                output1 << seq << endl;
                output1 << sym << endl;
                output1 << qual << endl;

                string line11 = temp1 + seq_bc + "-1";
                output2 << line11 << endl;
                output2 << seq1 << endl;
                output2 << sym1 << endl;
                output2 << qual1 << endl;
            }
            else
            {
                string line1 = temp;
                output1 << line1 << endl;
                output1 << seq << endl;
                output1 << sym << endl;
                output1 << qual << endl;

                string line11 = temp1;
                output2 << line11 << endl;
                output2 << seq1 << endl;
                output2 << sym1 << endl;
                output2 << qual1 << endl;
            }
        }
        else
        {
            string line1 = temp;
            output1 << line1 << endl;
            output1 << seq << endl;
            output1 << sym << endl;
            output1 << qual << endl;
     
            string line11 = temp1;
            output2 << line11 << endl;
            output2 << seq1 << endl;
            output2 << sym1 << endl;
            output2 << qual1 << endl;
        }
    }
    
    sam_close(fp1);
    free(bamHdr1);
    free(aln1);
    
    sam_close(fp2);
    free(bamHdr2);
    free(aln2);

    free(bamHdr3);
    free(aln3);
    sam_close(fp3);
}
