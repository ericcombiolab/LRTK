#include <htslib/sam.h>
#include <vector>
#include <unordered_map>
#include <thread>
#include <iostream>
#include <fstream>
#include "cmdline.h"
#include "ThreadPool/ThreadPool.h"

struct alignedRead
{
    unsigned long long read_id = 0, barcode_id = 0;
    unsigned int start_pos = 0, end_pos = 0;
    int isize = 0;
    bool isread1 = false;
};

struct readBasic
{
    unsigned long long read_id = 0;
    int read1len = 0, read2len = 0;
    unsigned int isize = 0;
    bool isDup = false;
};

struct barcodeBasic
{
    unsigned long long barcode_id = 0;
    unsigned int fragemnts_size = 0;
};

struct longFragment
{
    std::string reference;
    unsigned long long barcode = 0;
    unsigned int start_pos = 0, end_pos = 0, length = 0;
    int read_pairs = 0;
};

std::vector<longFragment> construct_fragment(std::string lastContig, std::vector<alignedRead> reads, int dist, int minL, int minR)
{
    std::vector<longFragment> ret;
    std::unordered_map<unsigned long long, std::vector<alignedRead>> barcode2reads;
    for (auto &&read : reads)
        barcode2reads[read.barcode_id].emplace_back(read);
    for (auto &&it : barcode2reads)
    {
        unsigned int start_pos = 0, end_pos = 0;
        int flag = 0, read_pair = 0;
        for (int i = 0; i < it.second.size(); ++i)
        {
            // flag == 0: start to search long fragments
            if (flag == 0)
            {
                int j = i + 1;
                for (; j < it.second.size(); ++j)
                {
                    if (it.second.at(j).read_id == it.second.at(i).read_id && it.second.at(j).isize + it.second.at(i).isize == 0 && ((it.second.at(j).isread1 && !it.second.at(i).isread1) || (!it.second.at(j).isread1 && it.second.at(i).isread1)))
                        break;
                }
                if (j != it.second.size())
                {
                    start_pos = it.second.at(i).start_pos;
                    end_pos = std::max(it.second.at(i).end_pos, it.second.at(j).end_pos);
                    flag = 1;
                    read_pair = 1;
                }
                else
                {
                    continue;
                }
            }
            // flag == 1: continue to search long fragments
            else if (flag == 1)
            {
                int j = i + 1;
                for (; j < it.second.size(); ++j)
                {
                    if (it.second.at(j).read_id == it.second.at(i).read_id && it.second.at(j).isize + it.second.at(i).isize == 0 && ((it.second.at(j).isread1 && !it.second.at(i).isread1) || (!it.second.at(j).isread1 && it.second.at(i).isread1)))
                        break;
                }
                if (j == it.second.size())
                {
                    continue;
                }
                else if (it.second.at(j).start_pos <= end_pos + dist)
                {
                    end_pos = std::max(end_pos, it.second.at(j).end_pos);
                    read_pair += 1;
                }
                else
                {
                    unsigned int fragL = end_pos - start_pos;
                    if (fragL >= minL && read_pair >= minR)
                    {
                        longFragment fragment;
                        fragment.barcode = it.first;
                        fragment.start_pos = start_pos;
                        fragment.end_pos = end_pos;
                        fragment.length = fragL;
                        fragment.reference = lastContig;
                        fragment.read_pairs = read_pair;
                        ret.emplace_back(fragment);
                    }
                    start_pos = it.second.at(i).start_pos;
                    end_pos = std::max(it.second.at(i).end_pos, it.second.at(j).end_pos);
                    read_pair = 1;
                }
            }
        }
        if (flag == 1)
        {
            unsigned int fragL = end_pos - start_pos;
            if (fragL >= minL && read_pair >= minR)
            {
                longFragment fragment;
                fragment.barcode = it.first;
                fragment.start_pos = start_pos;
                fragment.end_pos = end_pos;
                fragment.length = fragL;
                fragment.reference = lastContig;
                fragment.read_pairs = read_pair;
                ret.emplace_back(fragment);
            }
        }
    }
    reads.clear();
    barcode2reads.clear();
    return ret;
}

int main(int argc, char *argv[])
{
    // parse arguments
    cmdline::parser argParser;
    argParser.add<std::string>("bam", 'b', "alignment sorted bam file to reference genomes", true, "");
    argParser.add<double>("mapping-quality", 'q', "the threshold of mapping quality", false, 10);
    argParser.add<double>("identity", 'i', "the threshold of identity", false, 0.95);
    argParser.add<int>("distance", 'd', "the expanding distance of read pairs", false, 1000);
    argParser.add<int>("min-length", 'l', "the minimum length for long fragments", false, 1000);
    argParser.add<int>("min-read-pairs", 'r', "the minimum number of read pairs for long fragments", false, 2);
    argParser.add<int>("max-insert-size", 'm', "insert size cutoff", false, 1000);
    argParser.add<int>("threads", 't', "number of threads", false, 16);
    argParser.add<std::string>("output", 'o', "output path to the long fragments", true, "");
    argParser.parse_check(argc, argv);
    std::string bam = argParser.get<std::string>("bam");
    double mq = argParser.get<double>("mapping-quality");
    double idt = argParser.get<double>("identity");
    int dist = argParser.get<int>("distance");
    int minL = argParser.get<int>("min-length");
    int minR = argParser.get<int>("min-read-pairs");
    int insertS = argParser.get<int>("max-insert-size");
    int thread = argParser.get<int>("threads");
    std::string output = argParser.get<std::string>("output");

    // set necessary containers
    std::vector<alignedRead> reads;
    std::unordered_map<std::string, barcodeBasic> barcodes2basic;
    std::vector<std::string> barcodes;
    std::unordered_map<std::string, readBasic> reads2basic;
    std::unordered_map<std::string, unsigned int> reference2length;

    // set necessary variable for statistics
    double totalReadLength = 0, totalReferenceLength = 0, totalFragmentLength = 0, totalInsertSize = 0, NxFragment = 0;
    unsigned long cntFragments = 0;
    unsigned long long cntDup = 0, cntIS = 0, cntBC_constructed = 0;
    int N25Fragment = 0, N50Fragment = 0, N75Fragment = 0;
    bool N25flag = false, N50flag = false, N75flag = false;
    std::vector<unsigned int> fragmentsLength;
    std::vector<unsigned long long> bins(10);
    std::fstream fout(output, std::ios::out);
    fout << "ID\tBarcode\tReference_Name\tStart_Pos\tEnd_Pos\tLength\tRead_Pairs" << std::endl;

    // initialize multi-threading
    ThreadPool pool(thread);
    std::vector<std::future<std::vector<longFragment>>> results;

    const int supplementaryAlignment = 0x00000800, secondaryAlignment = 0x00000100;
    const int Read1 = 0x00000040, Read2 = 0x00000080, Dup = 0x00000400;
    samFile *bamFile = hts_open(bam.c_str(), "r");
    bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
    bam1_t *bamAln = bam_init1();
    std::string lastContig;
    unsigned long long cntLine = 0, cntRead = 0, cntBarcode = 0;
    std::cout << "Start parsing bam file ..." << std::endl;
    while (sam_read1(bamFile, bamHeader, bamAln) >= 0)
    {
        // set flags
        int flag = bamAln->core.flag;
        bool isSupplementary = flag & supplementaryAlignment, isSecondary = flag & secondaryAlignment;
        bool isRead1 = flag & Read1, isRead2 = flag & Read2, isDup = flag & Dup;
        // get read name
        std::string readName;
        if (bam_get_qname(bamAln) != NULL)
            readName = bam_get_qname(bamAln);
        if (readName.empty())
            continue;
        // update containers
        alignedRead read;
        if (reads2basic.find(readName) != reads2basic.end())
            read.read_id = reads2basic[readName].read_id;
        else
        {
            read.read_id = cntRead;
            reads2basic[readName].read_id = cntRead++;
        }
        if (isRead1)
            reads2basic[readName].read1len = bamAln->core.l_qseq;
        else
            reads2basic[readName].read2len = bamAln->core.l_qseq;
        if (isDup)
            reads2basic[readName].isDup = true;
        // get contigname
        std::string contigName;
        if (sam_hdr_tid2name(bamHeader, bamAln->core.tid) != NULL)
            contigName = bamHeader->target_name[bamAln->core.tid];
        if (contigName.empty())
            continue;
        reference2length[contigName] = bamHeader->target_len[bamAln->core.tid];
        // skip pcr duplicate, ambiguous reads, supplementary and secondary
        if (isDup || isSupplementary || isSecondary || (isRead1 && isRead2) || (!isRead1 && !isRead2))
            continue;
        // skip low mapping quality
        int mappingQuality = bamAln->core.qual;
        if (mappingQuality < mq)
            continue;
        // get alignment coordinates
        unsigned int pos = bamAln->core.pos + 1, endPos = bam_endpos(bamAln) + 1;
        if (endPos == pos + 1)
            continue;
        // get alignment length
        int alignedBases = 0, NM = 0;
        uint32_t *cigarPointer = bam_get_cigar(bamAln);
        for (int i = 0; i < bamAln->core.n_cigar; ++i)
        {
            char cigarOperator = bam_cigar_opchr(cigarPointer[i]);
            if (cigarOperator == 'M' || cigarOperator == 'I' || cigarOperator == 'D')
                alignedBases += bam_cigar_oplen(cigarPointer[i]);
        }
        if (alignedBases == 0)
            continue;
        // calculate blast identity
        uint8_t *tmpNM = bam_aux_get(bamAln, "NM");
        if (tmpNM != NULL)
            NM = bam_aux2i(tmpNM);
        double blastIdentity = 1.0 * (alignedBases - NM) / alignedBases;
        if (blastIdentity < idt)
            continue;
        // for insert size
        read.isize = bamAln->core.isize;
        if (read.isize >= 0)
            reads2basic[readName].isize = read.isize;
        else
            reads2basic[readName].isize = -read.isize;
        // get barcodes
        std::string barcode;
        uint8_t *tmpBC = bam_aux_get(bamAln, "BX");
        if (tmpBC != NULL)
            barcode = bam_aux2Z(tmpBC);
        // start parallel
        if (!lastContig.empty() && contigName.compare(lastContig))
        {
            results.emplace_back(pool.enqueue(construct_fragment, lastContig, reads, dist, minL, minR));
            reads.clear();
            if (results.size() >= 5)
            {
                for (auto &&result : results)
                {
                    std::vector<longFragment> fragments = std::move(result.get());
                    for (auto &&f : fragments)
                    {
                        cntFragments += 1;
                        totalFragmentLength += f.length;
                        fragmentsLength.emplace_back(f.length);
                        barcodes2basic[barcodes.at(f.barcode)].fragemnts_size += 1;
                        fout << cntFragments << "\t" << barcodes.at(f.barcode) << "\t" << f.reference << "\t" << f.start_pos << "\t" << f.end_pos << "\t" << f.length << "\t" << f.read_pairs << std::endl;
                    }
                }
                results.clear();
            }
        }
        if (contigName.compare(lastContig))
            lastContig = contigName;

        if (!barcode.empty())
        {
            read.start_pos = pos;
            read.end_pos = endPos;
            if (barcodes2basic.find(barcode) != barcodes2basic.end())
                read.barcode_id = barcodes2basic[barcode].barcode_id;
            else
            {
                read.barcode_id = cntBarcode;
                barcodes.emplace_back(barcode);
                barcodes2basic[barcode].barcode_id = cntBarcode++;
            }
            if (isRead1)
                read.isread1 = true;
            if (read.isize && read.isize <= insertS && read.isize >= -insertS)
                reads.emplace_back(read);
        }
        if (++cntLine % 10000000 == 0)
            std::cout << "\rParsed " << cntLine << " alignments ..." << std::flush;
    }
    results.emplace_back(pool.enqueue(construct_fragment, lastContig, reads, dist, minL, minR));
    reads.clear();
    for (auto &&result : results)
    {
        std::vector<longFragment> fragments = std::move(result.get());
        for (auto &&f : fragments)
        {
            cntFragments += 1;
            totalFragmentLength += f.length;
            fragmentsLength.emplace_back(f.length);
            barcodes2basic[barcodes.at(f.barcode)].fragemnts_size += 1;
            fout << cntFragments << "\t" << barcodes.at(f.barcode) << "\t" << f.reference << "\t" << f.start_pos << "\t" << f.end_pos << "\t" << f.length << "\t" << f.read_pairs << std::endl;
        }
    }
    results.clear();
    if (cntLine % 10000000)
        std::cout << "\rParsed " << cntLine << " alignments ...";
    std::cout << std::endl;
    bam_hdr_destroy(bamHeader);
    bam_destroy1(bamAln);
    sam_close(bamFile);
    fout.close();

    std::vector<unsigned int> inserts_vector;
    std::fstream foutis(output + ".insert_size", std::ios::out);
    foutis << "Insert_Size" << std::endl;
    for (auto &&read : reads2basic)
    {
        if (read.second.isDup)
            cntDup += 1;
        totalReadLength += (read.second.read1len + read.second.read2len);
        if (read.second.isize)
        {
            cntIS += 1;
            totalInsertSize += read.second.isize;
            inserts_vector.emplace_back(read.second.isize);
            foutis << read.second.isize << std::endl;
        }
    }
    std::sort(inserts_vector.begin(), inserts_vector.end(), std::greater<unsigned int>());
    foutis.close();
    for (auto &&ref : reference2length)
        totalReferenceLength += ref.second;
    for (auto &&bc : barcodes2basic)
    {
        if (bc.second.fragemnts_size)
        {
            cntBC_constructed += 1;
            if (bc.second.fragemnts_size <= 10)
                bins[bc.second.fragemnts_size - 1]++;
        }
    }
    std::sort(fragmentsLength.begin(), fragmentsLength.end(), std::greater<unsigned int>());
    for (auto &&f : fragmentsLength)
    {
        NxFragment += f;
        if (!N25flag && NxFragment > totalFragmentLength * 0.25)
        {
            N25flag = true;
            N25Fragment = f;
        }
        if (!N50flag && NxFragment > totalFragmentLength * 0.5)
        {
            N50flag = true;
            N50Fragment = f;
        }
        if (!N75flag && NxFragment > totalFragmentLength * 0.75)
        {
            N75flag = true;
            N75Fragment = f;
        }
    }
    std::cout << "\nTotal number of read pairs: " << cntRead << std::endl;
    std::cout << "Number of read pairs with non-zero insert size: " << cntIS << std::endl;
    std::cout << "Mean insert size: " << totalInsertSize / cntIS << std::endl;
    std::cout << "Median insert size: " << inserts_vector[inserts_vector.size()/2] << std::endl;
    std::cout << "Duplication ratio: " << cntDup << '/' << cntRead << '=' << 1.0 * cntDup / cntRead << std::endl;
    std::cout << "Total number of barcodes: " << cntBarcode << std::endl;
    std::cout << "Number of barcodes with fragemnts: " << cntBC_constructed << std::endl;
    std::cout << "Total number of fragments: " << cntFragments << std::endl;
    std::cout << std::endl;
    std::cout << "Total read length: " << totalReadLength << std::endl;
    std::cout << "Total fragment length: " << totalFragmentLength << std::endl;
    std::cout << "Total reference length: " << totalReferenceLength << std::endl;
    std::cout << std::endl;
    std::cout << "Read coverage (C): " << totalReadLength / totalReferenceLength << std::endl;
    std::cout << "Mean fragment length (U_FL): " << totalFragmentLength / cntFragments << std::endl;
    std::cout << "Fragment length max, N25, N50 (WU_FL), N75, min: " << fragmentsLength.at(0) << ", " << N25Fragment << ", " << N50Fragment << ", " << N75Fragment << ", " << fragmentsLength.at(fragmentsLength.size() - 1) << std::endl;
    std::cout << "Fragment coverage (C_F): " << totalFragmentLength / totalReferenceLength << std::endl;
    std::cout << "Read coverage per fragment (C_R): " << totalReadLength / totalFragmentLength << std::endl;
    std::cout << "Fragments per barcode (N_F/P): " << 1.0 * cntFragments / cntBC_constructed << std::endl;
    std::cout << "\nDistribution of fragments per barcode (<=10): \n";
    for (int i = 0; i < bins.size(); ++i)
        std::cout << i + 1 << "\t" << bins[i] << "\t" << 1.0 * bins[i] / cntBC_constructed << std::endl;
    return 0;
}