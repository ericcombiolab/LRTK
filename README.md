# LRTK: a unified ToolKit for Linked-Read sequencing data
# Overview
 LRTK consist of 4 sections: format conversion of FASTQ files, barcode-aware read alignment, quality control and downstream analysis modules. It can accept multiple linked-read format, including sequencing data from 10x genomics, stLFR, and TELL-Seq technologies. In the downstream analysis section, LRTK includes a variety of tools to the detect and phase multiple variations, such as SNPs, INDELs, and SVs.
 
![Main](https://user-images.githubusercontent.com/3699571/163749053-543ed7df-fb8f-4aa7-8c7e-3615b756e759.gif)
# Install through Bioconda
## bioconda install
(Please ensure channels are properly setup for bioconda before installing)
```
conda install lrtk (waiting)
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
If you would like to try the optional tools, such as LinkedSV, Valor, WhatSHAP, please make sure the program executables has been in one of the directories listed in the PATH environment variable (".bashrc").


In addition, the reference GRCH38 genome and related files can be found on the website: 
# Running examples
## module1: format conversion across diverse linked-read data format
```
LRTK FQCONVER -I1 NA12878_1.fastq -I2 NA12878_2.fastq -IT stLFR -O1 outFQ1 -O2 outFQ2 -OT ULRF -B BLstLFR -T 4 
```
## module2: unified barcode-aware alignment
```
LRTK ALIGN -BQ1 barcodedstLFRFQ1 -BQ2 barcodedstLFRFQ2 -FQ1 nobarcodedstLFRFQ1 -FQ2 nobarcodedstLFRFQ2 -R GRCH38.fa -O NA12878.bam -RG "@RG\tID:NA12878\tSM:NA12878" -P stLFR -T 4
```
## module3: variation calling and phasing
```
LRTK SNV -B NA12878.bam -R GRCH38 -A SAMTOOLS -T 4 -O NA12878.vcf
```
## module4: automatic pipeline to process multiple samples
```
LRTK WGS -SI SAMPLE_INFO -OD OUTDIR -DB DATABASE -T 4
```
# Troubleshooting
Please submit issues on the github page for LRTK.
