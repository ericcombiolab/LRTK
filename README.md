# LRTK: a unified ToolKit for Linked-Read sequencing data
# Overview
 LRTK is an easy-to-use tool kit to handle different linked-read sequencing data. It consist of multiple sections: format conversion of FASTQ files, barcode-aware read alignment, quality control and downstream analysis modules (SNP/INDEL/SV calling and phasing). It can accept multiple linked-read format, including sequencing data from 10x genomics, stLFR, and TELL-Seq technologies. 
 
![Main](https://user-images.githubusercontent.com/3699571/163749053-543ed7df-fb8f-4aa7-8c7e-3615b756e759.gif)
# Install through Bioconda
## bioconda install
(Please ensure channels are properly setup for bioconda before installing)
```
conda install lrtk
LRTK -h
```
## dependencies
LRTK utilizes python3 to execute the programs and requires some python packages, including numpy, pysam, scipy and sortedcontainers. some commonly used bioinformatics tools, such as SAMtools and minimap2, should be included in the path environment variable. Or you could directly install them through binconda.
```
conda install -c bioconda aquila
conda install -c bioconda bcftools
conda install -c bioconda bwa
conda install -c bioconda fastp
conda install -c bioconda freebayes
conda install -c bioconda hapcut2
conda install -c bioconda samtools
```
Some optional tools, such as LinkedSV, Valor and SpecHap, need to be installed from source codes and added into the PATH environment variable.

## example data in Zenodo
Please download the default database from Zenodo.
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
