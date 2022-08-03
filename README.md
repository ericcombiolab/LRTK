# LRTK: A unified and versatile toolkit for ana-lyzing linked-read sequencing data
# Overview
 LRTK is an easy-to-use tool kit to handle different linked-read sequencing data. It consist of multiple sections: format conversion of FASTQ files, barcode-aware read alignment, quality control and downstream analysis modules (SNP/INDEL/SV calling and phasing). It can accept multiple linked-read format, including sequencing data from 10x genomics, stLFR, and TELL-Seq technologies. 
 
![Main](https://user-images.githubusercontent.com/3699571/163749053-543ed7df-fb8f-4aa7-8c7e-3615b756e759.gif)
# Install through Bioconda
## bioconda install
(Please ensure channels are properly setup for bioconda before installing)
```
conda install -c bioconda lrtk
lrtk -h
```
## dependencies
LRTK is mainly implemented by python3 and requires python packages, such as numpy, pysam, scipy and sortedcontainers. The dependent tools will be installed automatically, including Aquila (Zhou et al., 2021), bcftools (Danecek et al., 2021), BWA (Li and Durbin, 2009), fastp (Chen et al., 2018), FreeBayes (Garrison and Marth, 2012), HapCUT2 (Edge et al., 2017) and SAMtools (Li et al., 2009). 
```
conda install -c bioconda aquila
conda install -c bioconda bcftools
conda install -c bioconda bwa
conda install -c bioconda fastp
conda install -c bioconda freebayes
conda install -c bioconda hapcut2
conda install -c bioconda samtools
```
Furthermore, LinkedSV (Fang et al., 2019), SpecHap (Yu et al., 2021) and VALOR2 (Karaoǧlanoǧlu et al., 2020) needed to be installed by the users because they are not supported by conda. 

## example data in Zenodo
The default database can be downloaded from Zenodo (https://zenodo.org/record/6792169).
The database directory should be organized as:
```
database
|-GRCH38
|-sonic
|-WhiteList
|-Uniqness_map
```
In addition, we provide multiple example linked-read sequencing data to run LRTK. We have included the small data sets (FQs) in github and stored the large data sets on Zenodo to generate the HTML report. 
```
example
|-FQs/
    |-Example.10x.R1.fq
    |-Example.10x.R2.fq
    |-Example.stLFR.R1.fq
    |-Example.stLFR.R2.fq
    |-Example.TELLSeq.index.fq
    |-Example.TELLSeq.R1.fq
    |-Example.TELLSeq.R2.fq
|-LargeFQs
    |-Example.10x.R1.fq
    |-Example.10x.R2.fq
    |-Example.stLFR.R1.fq
    |-Example.stLFR.R2.fq
    |-Example.TELLSeq.index.fq
    |-Example.TELLSeq.R1.fq
    |-Example.TELLSeq.R2.fq
```

# Running examples
To simplify the following description, we will define some commonly used variables here.
```
LRTK=path_to_lrtk
curP=path_to_working directory
DATABASE=path_to_default database

###input files
raw10xFQ1=${curP}"/FQs/Example.10x.R1.fq"
raw10xFQ2=${curP}"/FQs/Example.10x.R2.fq"
rawstLFRFQ1=${curP}"/FQs/Example.stLFR.R1.fq"
rawstLFRFQ2=${curP}"/FQs/Example.stLFR.R2.fq"
rawTELLSeqFQ1=${curP}"/FQs/Example.TELLSeq.R1.fq"
rawTELLSeqFQ2=${curP}"/FQs/Example.TELLSeq.R2.fq"
rawTELLSeqFQi=${curP}"/FQs/Example.TELLSeq.index.fq"

#output files
outFQ1=${curP}"/test/Example.R1.fq"
outFQ2=${curP}"/test/Example.R2.fq"
outBAM=${curP}"/test/Example.bam"
outVCF1=${curP}"/test/Example.small.variants.vcf"
outVCF2=${curP}"/test/Example.large.variants.vcf"
outVCF3=${curP}"/test/Example.small.variants.phased.vcf"

###database
BL10x=${DATABASE}"/WhiteList/white_list_10x_barcode.fa"
BLstLFR=${DATABASE}"/WhiteList/white_list_stlfr_barcode.fa"
BLTELLSeq=${DATABASE}"/WhiteList/4M-with-alts-february-2016.TELLSeq.txt"
GRCH38=${DATABASE}"/GRCH38/genome.fa"
SONIC=${DATABASE}"/sonic/GRCh38.sonic"
UNIQNESS=${DATABASE}"/Uniqness_map/"

```

## function 1: format conversion across diverse linked-read data formats
```
LRTK FQCONVER -I1 NA12878_1.fastq -I2 NA12878_2.fastq -IT stLFR -O1 outFQ1 -O2 outFQ2 -OT ULRF -B BLstLFR -T 4 
```
## function 2: unified barcode-aware alignment
```
LRTK ALIGN -BQ1 barcodedstLFRFQ1 -BQ2 barcodedstLFRFQ2 -FQ1 nobarcodedstLFRFQ1 -FQ2 nobarcodedstLFRFQ2 -R GRCH38.fa -O NA12878.bam -RG "@RG\tID:NA12878\tSM:NA12878" -P stLFR -T 4
```
## function 3: small variation calling
```
LRTK SNV -B NA12878.bam -R GRCH38 -A "SAMTOOLS" -T 4 -O NA12878.SNV.vcf
```
## function 4: large variation calling
```
LRTK SV -B NA12878.bam -R GRCH38 -A "Aquila" -T 4 -O NA12878.SV.vcf -V NA12878.SNV.raw.vcf -U path_to_UNIQNESS_database
```
## function 5: variation phasing
```
LRTK PHASE -B NA12878.bam -R GRCH38 -A "HapCUT2" -T 12 -V NA12878.SNV.vcf -O NA12878.SNV.phased.vcf
```
## function 6: automatic pipeline to process multiple samples
```
LRTK WGS -SI path_to_sample_info -OD path_to_outdir -DB path_to_database -RG "@RG\tID:NA12878\tSM:NA12878" -T 32
```
# Output files:

# Troubleshooting
If any question, please feel free to contact with me or my supervisor Dr.Eric Zhang.
(Chao Yang, email: cschaoyang@comp.hkbu.edu.hk)
(Lu Zhang, email: ericluzhang@hkbu.edu.hk)
