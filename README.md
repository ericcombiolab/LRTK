# LRTK: A unified and versatile toolkit for analyzing linked-read sequencing data
# Overview
 LRTK is an easy-to-use tool kit to handle linked-read sequencing data a unfrom 10x genomics, stLFR, and TELL-Seq technologies. It contains a suite of utilities to perform data simulation, format conversion, data preprocessing, barcode-aware read alignment, quality control, variant detection and phasing (SNP/INDEL/SV). In particular, LRTK is open-source and can gen-erate a HTML report to calculate the key parameters for library prepara-tion and summarize the quality statistics of sequenced reads.
 
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
In addition, we provide multiple example linked-read sequencing data to run LRTK. We have included the small data sets (FQs) in github and stored the large data sets on Zenodo  (https://zenodo.org/record/6792169) to generate the HTML report. 
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
###working environment
LRTK=path_to_lrtk
curP=path_to_working directory
DATABASE=path_to_default database
OUTDIR=${curP}"/WGS/"

###input files
raw10xFQ1=${curP}"/FQs/Example.10x.R1.fq"
raw10xFQ2=${curP}"/FQs/Example.10x.R2.fq"
rawstLFRFQ1=${curP}"/FQs/Example.stLFR.R1.fq"
rawstLFRFQ2=${curP}"/FQs/Example.stLFR.R2.fq"
rawTELLSeqFQ1=${curP}"/FQs/Example.TELLSeq.R1.fq"
rawTELLSeqFQ2=${curP}"/FQs/Example.TELLSeq.R2.fq"
rawTELLSeqFQi=${curP}"/FQs/Example.TELLSeq.index.fq"

barcoded10xFQ1=${curP}"/test/Example.R1.fq.sort.wb.fq"
barcoded10xFQ2=${curP}"/test/Example.R2.fq.sort.wb.fq"
nobarcoded10xFQ1=${curP}"/test/Example.R1.fq.sort.wob.fq"
nobarcoded10xFQ2=${curP}"/test/Example.R2.fq.sort.wob.fq"

Sinfo=${curP}"/sample.info"

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
## Commands for raw read and variant analysis
### function 1: example for format conversion across diverse linked-read data formats
```
$LRTK FQCONVER -I1 $raw10xFQ1 -I2 $raw10xFQ2 -IT 10x -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BL10x -T 4 
```
### function 2: example for unified barcode-aware alignment
```
$LRTK ALIGN -BQ1 $barcoded10xFQ1 -BQ2 $barcoded10xFQ2 -FQ1 $nobarcoded10xFQ1 -FQ2 $nobarcoded10xFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:Example\tSM:Example" -P 10x -T 4
```
### function 3: example for small variation calling
```
$LRTK SNV -B $outBAM -R $GRCH38 -A "SAMTOOLS" -T 4 -O $outVCF1
```
### function 4: example for large variation calling
```
$LRTK SV -B $outBAM -R $GRCH38 -A "Aquila" -T 4 -O $outVCF2 -V $outVCF1 -U $UNIQNESS
```
### function 5: example for variation phasing
```
$LRTK PHASE -B $outBAM -R $GRCH38 -A "HapCUT2" -T 4 -V $outVCF1 -O $outVCF3
```
## Commands for automatic pipeline
LRTK provides an easy-to-use automatic pipeline to handle the linked-read sequencing data from single or multiple samples. The user may only prepare the linked-read sequencing ﬁles (FASTQ format) and adequate computational resources, LRTK will run the whole pipeline and generate the final report. 
### function 6：example for single sample analysis
We show a simple example to process a single samples using the automatic pipeline. 
```
$LRTK WGS -SI $Sinfo -OD $OUTDIR -DB $DATABASE -RG "@RG\tID:Example\tSM:Example" -T 32
```
The sample information file should be organized as:
```
#Barcode	FQ1	FQ2	INDEXFQ	Linked-read_tech
Example_10x	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.10x.R1.fq	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.10x.R2.fq	-	10x
```
### function 7: example for multiple samples analysis
We show a simple example to process multiple samples using the automatic pipeline, simultaneously. 
```
$LRTK WGS -SI $Sinfo -OD $OUTDIR -DB $DATABASE -RG "@RG\tID:Example\tSM:Example" -T 32
```
The sample information file should be organized as:
```
#Barcode	FQ1	FQ2	INDEXFQ	Linked-read_tech
Example_10x	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.10x.R1.fq	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.10x.R2.fq	-	10x
Example_stLFR	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.stLFR.R1.fq	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.stLFR.R2.fq	-	stLFR
Example_TELLSeq	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.TellSeq.R1.fq	/tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.TellSeq.R2.fq /tmp/local/cschaoyang/SOFTWARE/LRTK/LRTK1.2/lrtk/example/LargeFQs/Example.TellSeq.index.fq TELLSeq
```

# Troubleshooting
If any question, please feel free to contact with me or my supervisor Dr.Eric Zhang.
(Chao Yang, email: cschaoyang@comp.hkbu.edu.hk)
(Lu Zhang, email: ericluzhang@hkbu.edu.hk)
