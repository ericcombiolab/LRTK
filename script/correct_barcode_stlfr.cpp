#include <htslib/sam.h>
# include <iostream> 
#include <string> 
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <gzstream.h>
#include <sstream>
using namespace std;

int main(int arg, char* argvs[])
{
    //argvs[1] barcode rightl list bc.fa
    //argvs[2] aligned barcode1 bc2.aln.sam
    //argvs[3] aligned barcode2 bc2.aln.sam
    //argvs[4] aligned barcode3 bc3.aln.sam
    //argvs[5] 1st reads without barcode. read1.fq
    //argvs[6] 2nd reads withtou barcode read2.fq
    //argvs[7] number of mismatch allowed NM
    //argvs[8] output file name for the read 1
    //argvs[9] output file name for the read 2
    cout << "Start output reads" << endl;
    unordered_map <string, string> ref_bc;
    fstream barcode(argvs[1], ios::in);
    string temp;
    while (getline(barcode, temp))
    {
        string refName = temp.substr(1, temp.length());
        getline(barcode, temp);
        ref_bc.insert({ refName, temp });
    }
    barcode.close();

    samFile* fp1 = hts_open(argvs[2], "r");
    samFile* fp2 = hts_open(argvs[3], "r");
    samFile* fp3 = hts_open(argvs[4], "r");
    bam_hdr_t* bamHdr1 = sam_hdr_read(fp1);
    bam1_t* aln1 = bam_init1();
    bam_hdr_t* bamHdr2 = sam_hdr_read(fp2);
    bam1_t* aln2 = bam_init1();
    bam_hdr_t* bamHdr3 = sam_hdr_read(fp3);
    bam1_t* aln3 = bam_init1();
    int NM_patience = atoi(argvs[7]);

    string temp1, temp2;
    string seq1, sym1, qual1, seq2, sym2, qual2;
    igzstream read1(argvs[5]);
    igzstream read2(argvs[6]);
    fstream output1(argvs[8], ios::out);
    fstream output2(argvs[9], ios::out);
    const int unmapped = 0x4;
    unsigned long long cnt_unmapped = 0, cnt_filtered = 0;
    unordered_set<string> bcs_out;

    while (sam_read1(fp1, bamHdr1, aln1) == 0 && sam_read1(fp2, bamHdr2, aln2) == 0 && sam_read1(fp3, bamHdr3, aln3) == 0)
    {
        getline(read1, temp1);
        getline(read2, temp2);
        istringstream iss(temp1);
        string readName1;
        iss >> readName1;
        istringstream iss2(temp2);
        string readName2;
        iss2 >> readName2;
        readName1 = readName1.substr(0, readName1.find("/1"));
        readName2 = readName2.substr(0, readName2.find("/2"));

        getline(read1, seq1);
        getline(read1, sym1);
        getline(read1, qual1);
        getline(read2, seq2);
        getline(read2, sym2);
        getline(read2, qual2);
        seq2 = seq2.substr(0, 100);
        qual2 = qual2.substr(0, 100);

        int flag1 = aln1->core.flag;
        int flag2 = aln2->core.flag;
        int flag3 = aln3->core.flag;

        if ((flag1 & unmapped) || (flag2 & unmapped) || (flag3 & unmapped))
        {
            cnt_unmapped += 1;
            output1 << readName1 << '\n' << seq1 << '\n' << sym1 << '\n' << qual1 << '\n';
            output2 << readName2 << '\n' << seq2 << '\n' << sym2 << '\n' << qual2 << '\n';
            continue;
        }

        uint8_t* NM1 = bam_aux_get(aln1, "NM");
        int64_t nm1 = bam_aux2i(NM1);
        uint8_t* NM2 = bam_aux_get(aln2, "NM");
        int64_t nm2 = bam_aux2i(NM2);
        uint8_t* NM3 = bam_aux_get(aln3, "NM");
        int64_t nm3 = bam_aux2i(NM3);

        if (nm1 + nm2 + nm3 <= NM_patience)
        {
            std::string refname1 = bamHdr1->target_name[aln1->core.tid];
            std::string refname2 = bamHdr2->target_name[aln2->core.tid];
            std::string refname3 = bamHdr3->target_name[aln3->core.tid];
            auto bc_loc1 = ref_bc.find(refname1);
            auto bc_loc2 = ref_bc.find(refname2);
            auto bc_loc3 = ref_bc.find(refname3);
            string seq_bc = " BX:Z:" + bc_loc1->second + bc_loc2->second + bc_loc3->second + "-1";
            bcs_out.insert(bc_loc1->first + "_" + bc_loc2->first + "_" + bc_loc3->first);
            readName1 += seq_bc;
            readName2 += seq_bc;
        }
        else
            cnt_filtered += 1;
        output1 << readName1 << '\n' << seq1 << '\n' << sym1 << '\n' << qual1 << '\n';
        output2 << readName2 << '\n' << seq2 << '\n' << sym2 << '\n' << qual2 << '\n';
    }
    cout << cnt_unmapped << " reads unmapped, " << cnt_filtered << " reads filtered (mismatch > " << NM_patience << "), unique barcodes " << bcs_out.size() << endl;
    bam_hdr_destroy(bamHdr1);
    bam_hdr_destroy(bamHdr2);
    bam_hdr_destroy(bamHdr3);
    bam_destroy1(aln1);
    bam_destroy1(aln2);
    bam_destroy1(aln3);
    sam_close(fp1);
    sam_close(fp2);
    sam_close(fp3);
    output1.close();
    output2.close();

    return 0;
}
