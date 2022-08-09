# LRTK: A unified and versatile toolkit for analyzing linked-read sequencing data
# Overview
 LRTK is an easy-to-use toolkit to process linked-read sequencing data from 10x genomics, stLFR, and TELL-Seq technologies. It contains a suite of utilities to perform data simulation, format conversion, data preprocessing, barcode-aware read alignment, quality control, variant detection and phasing (SNV/INDEL/SV). In particular, LRTK is open-source and can generate a HTML report to calculate the key parameters for library preparation and summarize the quality statistics of sequenced reads.
 
![Main](https://user-images.githubusercontent.com/3699571/163749053-543ed7df-fb8f-4aa7-8c7e-3615b756e759.gif)
# Install using Bioconda
## bioconda install
(Please ensure channels are properly setup for bioconda before installing)
```
conda install -c bioconda lrtk
lrtk -h
```
## dependencies
LRTK is mainly implemented by python3 and requires several python packages, such as numpy, pysam, scipy and sortedcontainers. Some of them will be installed using conda automatically, including Aquila (Zhou et al., 2021), bcftools (Danecek et al., 2021), BWA (Li and Durbin, 2009), fastp (Chen et al., 2018), FreeBayes (Garrison and Marth, 2012), HapCUT2 (Edge et al., 2017) and SAMtools (Li et al., 2009). 
Furthermore, <a href="https://github.com/WGLab/LinkedSV"> LinkedSV </a> (Fang et al., 2019), <a href="https://github.com/deepomicslab/SpecHap"> SpecHap </a> (Yu et al., 2021) and <a href="https://github.com/BilkentCompGen/valor"> VALOR2 </a> (Karaoǧlanoǧlu et al., 2020) that are not supported by conda are needed to be installed by users. 



## Example data in Zenodo
The required database can be downloaded from Zenodo (https://zenodo.org/record/6792169).
The directory is organized as:
```
database
|-GRCH38
|-sonic
|-WhiteList
|-Uniqness_map
```
In addition, we provide several examples to test LRTK on different linked-read sequencing technologies. We have included a small dataset (FQs) and a large dataset (LargeFQs) on Zenodo (https://zenodo.org/record/6792169). The small dataset is used to quickly test the raw read analysis module.The large data set can be used to test all the fucntions and generate the final HTML report.
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
    |-Example.large.10x.R1.fq
    |-Example.large.10x.R2.fq
    |-Example.large.stLFR.R1.fq
    |-Example.large.stLFR.R2.fq
    |-Example.large.TELLSeq.index.fq
    |-Example.large.TELLSeq.R1.fq
    |-Example.large.TELLSeq.R2.fq
```

# Running examples
## Help information
Type "lrtk -h" and the help information will be printed.
```
lrtk -h 
```
```
usage: lrtk version 1.5

Linked Reads ToolKit

positional arguments:
  {checkENV,MKFQ,FQCONVER,ALIGN,SNV,SV,PHASE,WGS}
    checkENV            Check the working environment
    MKFQ                Simulate linked-reads
    FQCONVER            Convert FASTQ formats
    ALIGN               Align reads to the reference genome
    SNV                 Detect SNVs and INDELs
    SV                  Detect structural variations
    PHASE               Phase germline variations
    WGS                 Run the whole pipeline

optional arguments:
  -h, --help            show this help message and exit

```
To simplify the following descriptions in running examples, we will firstly define some commonly used variables here.
```
LRTK=`which lrtk`
curP=`pwd`
DATABASE=path_to_default database
OUTDIR="./WGS/"

###input files
#We have included a small data sets (FQs) and a large data sets (LargeFQs) on Zenodo (https://zenodo.org/record/6792169).
#The small data set can be used to quickly test the raw read analysis module.
#The large data set can be used to test all the commands and generate the final HTML report.

raw10xFQ1="./example/FQs/Example.10x.R1.fq"
raw10xFQ2="./example/FQs/Example.10x.R2.fq"
rawstLFRFQ1="./example/FQs/Example.stLFR.R1.fq"
rawstLFRFQ2="./example/FQs/Example.stLFR.R2.fq"
rawTELLSeqFQ1="./example/FQs/Example.TELLSeq.R1.fq"
rawTELLSeqFQ2="./example/FQs/Example.TELLSeq.R2.fq"
rawTELLSeqFQi="./example/FQs/Example.TELLSeq.index.fq"

raw10xFQ1_large="./example/LargeFQs/Example.large.10x.R1.fq"
raw10xFQ2_large="./example/LargeFQs/Example.large.10x.R2.fq"
rawstLFRFQ1_large="./example/LargeFQs/Example.large.stLFR.R1.fq"
rawstLFRFQ2_large="./example/LargeFQs/Example.large.stLFR.R2.fq"  
rawTELLSeqFQ1_large="./example/LargeFQs/Example.large.TELLSeq.R1.fq"  
rawTELLSeqFQ2_large="./example/LargeFQs/Example.large.TELLSeq.R2.fq"  
rawTELLSeqFQi_large="./example/LargeFQs/Example.large.TELLSeq.index.fq" 

### intermediate results
barcoded10xFQ1="./Example.R1.fq.sort.wb.fq"  
barcoded10xFQ2="./Example.R2.fq.sort.wb.fq"  
nobarcoded10xFQ1="./Example.R1.fq.sort.wob.fq"  
nobarcoded10xFQ2="./Example.R2.fq.sort.wob.fq"  

### should be prepared by users to test the WGS command
Sinfo="./sample.info"  

###output files
outFQ1="./Example.R1.fq"
outFQ2="./Example.R2.fq"
outBAM="./Example.bam"
outVCF1="./Example.small.variants.vcf"
outVCF2="./Example.large.variants.vcf"
outVCF3="./Example.small.variants.phased.vcf"

###database
BL10x="./database/WhiteList/white_list_10x_barcode.fa"
BLstLFR="./database/WhiteList/white_list_stlfr_barcode.fa"
BLTELLSeq="./database/WhiteList/4M-with-alts-february-2016.TELLSeq.txt"
GRCH38="./database/GRCH38/genome.fa"
SONIC="./database/sonic/GRCh38.sonic"
UNIQNESS="./database/Uniqness_map/"
```
## Commands for raw read and variant analysis

### function 1: example for linked-read simulation
```
$LRTK MKFQ -CF "./example/FQs/simulation/diploid_config" -IT 10x
```
*parameters:

-CF/--config_file: The path to config_files for linked-read simulation

-IT/--input_type: Platform. Users can choose from 10x or stLFR.

We have prepared Two examples of config file (config1.txt and config2.txt; config1.txt is illustrated here) in the diploid_config folder.
```
# path of template fasta file
Path_Fastahap1=/path_to/example/FQs/simulation/resource/hap1.fa
Path_Fastahap2=/path_to/example/FQs/simulation/resource/hap2.fa
# number of processors for parallele
processors=50
# coverage for long fragment
CF=15
# coverage for short reads
CR=0.2
# the average number of molecules for each droplet
N_FP=16
# the average length for long fragment (Kb)
Mu_F=20
# length of short reads (bp)
SR=100
# fast mode ('Y' or 'N'; only simulate uniform sequencing quality)
Fast_mode=N
# simulate sequencing error ('Y') or not ('N')
Seq_error=Y
# sequencing error rate
Error_rate=0.01
#path to sequencing error profile
Path_Seq_qual=/path_to/example/FQs/simulation/resource/error_profile_reads.txt
#path to barcode error profile
Path_Barcode_qual=/path_to/example/FQs/simulation/resource/error_profile_barcode.txt
# mean of insert size for short reads (bp)
Mu_IS=400
# standard deviation of insert size for short reads (bp)
Std_IS=10
# barcode list
Path_barcodepool=/path_to/example/FQs/simulation/resource/barcodelist.txt
# Haploid (Hap=1) or Diploid (Hap=2)
Hap=2
Barcode_Length=54
```
line2 and line3: Path_Fastahap1 and Path_Fastahap2, the two haploid reference serquences. 

line5: processors, the maximum number of CPUs allowed

line7: CF, coverage of long DNA fragments

line9: CR, covergae of short reads for each fragment

line11: N_FP, average number of fragments for each droplet

line13: Mu_F, average length for long DNA fragment (Kb)

line15: SR, length of short reads (bp)

line21: Error_rate, sequencing error rate

line27: Mu_IS, the average of insert size for short reads (bp)

line29: Std_IS, standard deviation of insert size for short reads (bp)

line33: Hap, Haploid (Hap=1) or Diploid (Hap=2)

### function 2: example for format conversion across different linked-read data formats
```
$LRTK FQCONVER -I1 $raw10xFQ1 -I2 $raw10xFQ2 -IT 10x -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BL10x -T 4 
```
*parameters:

-I1/--input_fastq1: Input fastq file (uncompressed FASTQ format) for the first read of paired-end linked-read sequencing data.

-I2/--input_fastq2: Input fastq file (uncompressed FASTQ format) for the second read of paired-end linked-read sequencing data.

-ID/--index_fastq: Input index file (uncompressed FASTQ format) for paired linked-read sequencing data.

-IT/--input_type: Platforms. Users can choose from 10x,stLFR or TELLSeq.

-O1/--output_fastq1: Output fastq file for the first read of paired-end linked-read sequencing data.

-O2/--output_fastq2: Output fastq file for the second read of paired-end linked-read sequencing data.

-OT/--output_type: Output fastq format. The unified linked read format (ULRF) is recommended.

-B/--barcodes: Path to barcode whitlist.

-F/--filter: Yes/No. "Yes" indicates that LRTK will use fastp to filter reads.

-S/--sort: Yes/No. "Yes" indicates that LRTK will sort the reads based on barcodes.

-T/--threads: default = 1, the number of threads used for BWA and SAMtools.

### function 3: example for unified barcode-aware alignment
```
$LRTK ALIGN -BQ1 $barcoded10xFQ1 -BQ2 $barcoded10xFQ2 -FQ1 $nobarcoded10xFQ1 -FQ2 $nobarcoded10xFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:Example\tSM:Example" -P 10x -T 4
```
*parameters:

-BQ1/--input_barcoded_fastq1: Input fastq file (uncompressed FASTQ format) for the first read of paired linked-read sequencing data (with barcode).

-BQ2/--input_barcoded_fastq2: Input fastq file (uncompressed FASTQ format) for the second read of paired linked-read sequencing data (with barcode).

-FQ1/--input_fastq1:  Input fastq file (uncompressed FASTQ format) for the first read of paired linked-read sequencing data (without barcode).

-FQ2/--input_fastq2:  Input fastq file (uncompressed FASTQ format) for the second read of paired linked-read sequencing data (without barcode).

-RG/--read_group:  Full read group string (e.g. '@RG\tID:foo\tSM:bar')

-R/--reference:  The recommended "GRCH38/genome.fa" is the reference fasta file downloaded from Zenodo.

-O/--outfile:  The output alignment file.  

-S/--sort: Users can choose from (Yes, No). "Yes" indicates that LRTK will use samtools to sort alignment files based on genomic coordinate.

-M/--mark_duplication: Users can choose from (Yes, No). "Yes" indicates that LRTK will use picard to mark the duplicated reads using barcode information.

-P/--platform: Input sequencing technology. Users can choose from (10x,stLFR,TELLSeq).

-T/--threads: default = 1, this determines the number of threads used for ema, bwa and samtools.

### function 4: example for small variation calling
```
$LRTK SNV -B $outBAM -R $GRCH38 -A "FreeBayes" -T 4 -O $outVCF1
```
*parameters

-B/--bam: The alignment file (.bam) obtained from ALIGN function.  

-R/--reference: Reference genome.

-A/--application: The SNV/INDEL caller (FreeBayes, Samtools or GATK).

-T/--threads: The number of threads used in SNV/INDEL caller (default: 1).

-O/--outfile: Output VCF file.

### function 5: example for structural variation calling
```
$LRTK SV -B $outBAM -R $GRCH38 -A "Aquila" -T 4 -O $outVCF2 -V $outVCF1 -U $UNIQNESS
```
*parameters

-B/--bam: The alignment file (.bam) obtained from ALIGN function.  

-R/--reference: Reference genome.

-A/--application: The SV caller (Aquila, LinkedSV or VALOR).

-T/--threads: The number of threads used in SNV/INDEL caller (default: 1).

-O/--outfile: Output VCF file.

-U/--uniqness: "Uniqness_map/" is a required database for Aquila, which can be downloaded from Zenodo.

-S/--sonic: "sonic/" is a required database for VALOR, which can be downloaded from Zenodo.

-V/--vcf: VCF file generated from ```SNV``` function.


### function 6: example for variant phasing
```
$LRTK PHASE -B $outBAM -R $GRCH38 -A "HapCUT2" -T 4 -V $outVCF1 -O $outVCF3
```
*parameters

-B/--bam: The alignment file (.bam) obtained from ALIGN function.  

-R/--reference: Reference genome.

-A/--application: The variant phasing tool (HapCUT2, WhatsHap or SpecHap).

-T/--threads: The number of threads used for variant phasing tool (default: 1). 

-O/--outfile: Output phased VCF file.

-V/--vcf: VCF with variants to phase.

## Commands for automatic pipeline
LRTK provides an easy-to-use automatic pipeline to process linked-read sequencing from single or multiple samples. Users only need to prepare FASTQ files, LRTK will run the whole pipeline and generate the final report. 

### function 7：example for single sample analysis
We show a simple example to process a single samples using the automatic pipeline. 
```
$LRTK WGS -SI $Sinfo -OD $OUTDIR -DB $DATABASE -RG "@RG\tID:Example\tSM:Example" -T 32
```
*parameters
-OD/--outdir: The output directory.

-DB/--database: The ```database``` contains reference genome and barcode whitelist, which can be downloaded from Zenodo.

-RG/--read_group: Read group (e.g. '@RG\tID:foo\tSM:bar').

-T/--threads: The number of threads used for variant phasing tool (default: 1).

-SI/--sample_info: The path to Sinfo (sample information file).

The Sinfo (tab-separated) should be prepared as: 
```
#Barcode	FQ1	FQ2	INDEXFQ	Linked-read_tech
Example_10x	/path_to/Example.large.10x.R1.fq	/path_to/Example.large.10x.R2.fq	-	10x
```
### function 8: example for multiple samples analysis
We show a simple example to process multiple samples using the automatic pipeline, simultaneously. 
```
$LRTK WGS -SI $Sinfo -OD $OUTDIR -DB $DATABASE -RG "@RG\tID:Example\tSM:Example" -T 32
```
*parameters
-OD/--outdir: The output directory.

-DB/--database: The ```database``` contains reference genome and barcode whitelist, which can be downloaded from Zenodo.

-RG/--read_group: Read group (e.g. '@RG\tID:foo\tSM:bar').

-T/--threads: The number of threads used for variant phasing tool (default: 1).

-SI/--sample_info: The path to Sinfo (sample information file).

The Sinfo (tab-separated) should be prepared as:
```
#Barcode	FQ1	FQ2	INDEXFQ	Linked-read_tech
Example_10x	/path_to/Example.large.10x.R1.fq	/path_to/Example.large.10x.R2.fq	-	10x
Example_stLFR	/path_to/Example.large.stLFR.R1.fq	/path_to/Example.large.stLFR.R2.fq	-	stLFR
Example_TELLSeq	/path_to/Example.large.TellSeq.R1.fq	/path_to/Example.large.TellSeq.R2.fq /path_to/Example.large.TellSeq.index.fq TELLSeq
```
# examples of report: 
lrtk creates reports in the HTML format.

<a href="https://github.com/ericcombiolab/LRTK/blob/main/demo_report/LRTK.report.html">demo report</a> 

# Troubleshooting
If any question, please feel free to contact with me or my supervisor Dr.Eric Zhang.
(Chao Yang, email: cschaoyang@comp.hkbu.edu.hk)
(Lu Zhang, email: ericluzhang@hkbu.edu.hk)
