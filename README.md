# LRTK: A platform agnostic toolkit for linked-read analysis of both human genomes and metagenomes
# Overview
Linked-Read ToolKit (LRTK), is a unified and versatile toolkit to process human and metagenomic linked-read sequencing data from different linked-read sequencing technologies, including 10x Genomics, single-tube long fragment read (stLFR) and transposase enzyme linked long-read sequencing (TELL-Seq). LRTK provides functions to perform read cloud assembly, barcode-aware read alignment, reconstruction of long DNA fragments, taxonomic classification and quantification, and genomic variant calling and phasing. LRTK also has the ability to perform automatically analyze from raw sequencing data to downstream analysis and support the analysis of multiple samples in parallel. In addition, LRTK could produce publication-ready visualizations and generates reproducible reports, summarizing the key parameters at multiple checkpoints such as library preparation. 
 example
![Main](https://github.com/CicyYeung/LinkedReadToolKit/blob/main/script/HTML/img/workflow.png)
# Install using Bioconda
## Bioconda install
(Please ensure channels are properly setup for bioconda before installing)
```
conda install -c bioconda lrtk
lrtk -h
```
## Dependencies
LRTK is mainly implemented by python3 and requires several python packages, such as numpy, pandas, pysam, scipy, sklrean, snakemake, torch, and sortedcontainers. Most of them could be installed using conda automatically, including Aquila (Zhou et al., 2021), bcftools (Danecek et al., 2021), BWA (Li and Durbin et al., 2009), fastp (Chen et al., 2018), FreeBayes (Garrison and Marth et al., 2012), HapCUT2 (Edge et al., 2017), inStrain (Olm M R et al., 2021), Pangaea (Zhang Z et al., 2022), SAMtools (Li et al., 2009). spades (<=3.15), jellyfish(2.3.0), WhatsHap (Patterson M et al., 2015). 
Furthermore, <a href="https://github.com/WGLab/LinkedSV"> LinkedSV </a> (Fang et al., 2019) and <a href="https://github.com/BilkentCompGen/valor"> VALOR2 </a> (Karaoǧlanoǧlu et al., 2020) that are not supported by conda are needed to be installed by users. 

We have tested LRTK with Aquila==1.0.0, bcftools==1.8, bwa==0.7.17, fastp==0.23.2, freebayes==0.9.21, gatk==3.8, hapcut2==1.3.3, samtools==1.6, VALOR==2.1.5, whatshap==1.2.1.   


## Database and demo data 
The required database can be downloaded from Google Drive (https://drive.google.com/drive/folders/1XPW2avL_LZAt5yIh9tb35jZ5GfCSj7eQ).
In addition, we provide several examples to test LRTK on different linked-read sequencing technologies. We have included a human genome dataset (FQs) and a metagenomic dataset on Google Drive (https://drive.google.com/drive/folders/1XPW2avL_LZAt5yIh9tb35jZ5GfCSj7eQ).

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
  {MKFQ, FQCONVER, ALIGN, RLF, SNV, SV, PHASE, ASSEMBLY, WGS, MWGS}
    MKFQ                Simulate linked-reads
    FQCONVER            Convert linked FASTQ formats
    ALIGN               Barcode aware alignment to map reads to the reference genome
    RLF                 Reconstruct the long DNA fragment
    SNV                 Detect SNVs and INDELs
    SV                  Detect structural variations
    PHASE               Phase germline variations
    ASSEMBLY            Read-cloud metagenome assembly
    WGS                 Run the human genome sequencing analysis pipeline
    MWGS                Run the metagenome sequencing analysis pipeline

optional arguments:
  -h, --help            show this help message and exit
```

## Commands for raw read and variant analysis

### function 1: linked-read simulation
```
LRTK MKFQ -CF "/path_to/diploid_config" -IT stLFR
```
*parameters:

-CF/--config_file: The path to config_files for linked-read simulation

-IT/--input_type: Platform. Users can choose from 10x or stLFR.

We have prepared Two examples of config file (config1.txt and config2.txt; config1.txt is illustrated here) in the diploid_config folder.

### function 2: unified linked-read format conversion
```
LRTK FQCONVER -I1 /path_to/IN_FQ1 -I2 /path_to/INF_Q2 -IT 10x -O1 /path_to/OUT_FQ1 -O2 /path_to/OUT_FQ2 -B /path_to/BARCODE_WHITELIST -T 4 
```
*parameters:

-I1/--input_fastq1: Input fastq file (uncompressed FASTQ format) for the first read of paired-end linked-read sequencing data.

-I2/--input_fastq2: Input fastq file (uncompressed FASTQ format) for the second read of paired-end linked-read sequencing data.

-ID/--index_fastq: Input index file (uncompressed FASTQ format) for paired linked-read sequencing data. This files only contains information about barcode sequence.

-IT/--input_type: Platforms. Users can choose from 10x,stLFR or TELLSeq.

-O1/--output_fastq1: Output fastq file for the first read of paired-end linked-read sequencing data.

-O2/--output_fastq2: Output fastq file for the second read of paired-end linked-read sequencing data.

### function 3: barcode-aware alignment
```
LRTK ALIGN -FQ1 /path_to/IN_FQ1 -FQ2 /path_to/IN_FQ2 -R /path_to/REFERENCE -O /path_to/OUT_BAM -RG "@RG\tID:Example\tSM:Example" -P 10x -T 4
```
*parameters:

-FQ1/--input_fastq1:  Input fastq file (uncompressed FASTQ format) for the first read of paired linked-read sequencing data (without barcode).

-FQ2/--input_fastq2:  Input fastq file (uncompressed FASTQ format) for the second read of paired linked-read sequencing data (without barcode).

-RG/--read_group:  Full read group string (e.g. '@RG\tID:foo\tSM:bar')

-R/--reference:  The recommended "GRCH38/genome.fa" is the reference fasta file downloaded from Zenodo.

-O/--outfile:  The output alignment file.  

-P/--platform: Input sequencing technology. Users can choose from (10x,stLFR,TELLSeq).

### function 4: Reconstruction of long DNA fragment
```
LRTK RLF -B /path_to/IN_BAM -D 200000 -O /path_to/OUTFILE
```
*parameters:

-B/--bam: The alignment file (.bam) obtained using barcode aware approach.  

-D/--distance: the distance for seed extension 

-O/--outfile:  The reconstructed fragments.  

### function 5: small variation calling
```
LRTK SNV -B /path_to/IN_BAM -R /path_to/REFERENCE -A "FreeBayes" -T 4 -O /path_to/OUT_VCF
```
*parameters

-B/--bam: The alignment file (.bam) obtained from ALIGN function.  

-R/--reference: Reference genome.

-A/--application: The SNV/INDEL caller (FreeBayes, Samtools or GATK).

-O/--outfile: Output VCF file.

### function 6: structural variation calling
```
LRTK SV -B /path_to/IN_BAM -R /path_to/REFERENCE -A "Aquila" -T 4 -O /path_to/OUT_VCF -V /path_to/IN_VCF -U /path_to/DATABASE_UNIQNESS
```
*parameters

-B/--bam: The alignment file (.bam) obtained from ALIGN function.  

-R/--reference: Reference genome.

-A/--application: The SV caller (Aquila, LinkedSV or VALOR).

-O/--outfile: Output VCF file.

-U/--uniqness: "Uniqness_map/" is a required database for Aquila, which can be downloaded from Google Drive.

-S/--sonic: "sonic/" is a required database for VALOR, which can be downloaded from Google Drive.

-V/--vcf: VCF file generated from ```SNV``` function.


### function 7:  variant phasing
```
LRTK PHASE -B /path_to/IN_BAM -R /path_to/REFERENCE -A "HapCUT2" -V /path_to/IN_VCF -O /path_to/OUT_VCF
```
*parameters

-B/--bam: The alignment file (.bam) obtained from ALIGN function.  

-R/--reference: Reference genome.

-A/--application: The variant phasing tool (HapCUT2, WhatsHap or SpecHap).

-T/--threads: The number of threads used for variant phasing tool (default: 1). 

-O/--outfile: Output phased VCF file.

-V/--vcf: VCF with variants to phase.

### function 8:  metagenome assembly
```
LRTK ASSEMBLY -FQ1 /path_to/IN_FQ1 -FQ2 /path_to/IN_FQ2 -MS /path_to/METASPADES_CONTIG -AL /path_to/ATHENA_LOCAL_CONTIG -AH /path_to/ATHENA_HYBRID_CONTIG -LT LOW_ABD_CUT -O OUTFILE
```
*parameters

  -FQ1/--fq1:        Input FASTQ1 file (uncompressed FASTQ format).
  
  -FQ2/--fq2:        Input FASTQ2 file (uncompressed FASTQ format).
  
  -MS/--metaspades:   assembled contigs from metaspades.
  
  -AL/--athena_l:     local assembled contigs from athena.
  
  -AH/--athena_h:     hybrid assembled contigs from athena.
  
  -LT/--low_abd_cut:  coverage for low abundance contigs.
  
  -O/ --outfile:      the final assembled contigs.

## Commands for automatic analysis pipeline
LRTK provides an easy-to-use automatic pipeline to process linked-read sequencing from single or multiple samples. Users only need to prepare FASTQ files, LRTK will run the whole pipeline and generate the final report. We show simple examples to process human and metagenome sequencing data using the automatic pipeline. 

### function 9：human genome workflow
```
LRTK WGS -SI /path_to/SAMPLE_INFO -OD /path_to/OUTDIR -DB /path_to/DATABASE -RG "@RG\tID:Example\tSM:Example" 
```
*parameters
-OD/--outdir: The output directory.

-DB/--database: The ```database``` contains reference genome and barcode whitelist, which can be downloaded from Zenodo.

-RG/--read_group: Read group (e.g. '@RG\tID:foo\tSM:bar').

-SI/--sample_info: The path to Sinfo (sample information file).

The Sinfo (tab-separated) should be prepared as: 
```
#Barcode	FQ1	FQ2	INDEXFQ	Linked-read_tech
Example_10x	/path_to/Example.large.10x.R1.fq	/path_to/Example.large.10x.R2.fq	-	10x
```
### function 10: metagenome workflow
```
LRTK MWGS -SI /path_to/SingleSample_info -MI /path_to/MultipleSample_info -OD /path_to/OUTDIR -DB /path_to/DATABASE -RG "@RG\tID:foo\tSM:bar"
```
*parameters

-OD/--outdir: The output directory.

-DB/--database: The ```database``` contains reference genome and barcode whitelist, which can be downloaded from Google Drive.

-RG/--read_group: Read group (e.g. '@RG\tID:foo\tSM:bar').

-SI/--sample_info: The path to Sinfo (sample information file).

-MI/--multi_info: The path to Sinfo (sample information file).

The Sinfo (tab-separated) should be prepared as:
```
#Barcode	FQ1	FQ2	INDEXFQ	Linked-read_tech
Example_10x	/path_to/Example.large.10x.R1.fq	/path_to/Example.large.10x.R2.fq	-	10x
Example_stLFR	/path_to/Example.large.stLFR.R1.fq	/path_to/Example.large.stLFR.R2.fq	-	stLFR
Example_TELLSeq	/path_to/Example.large.TellSeq.R1.fq	/path_to/Example.large.TellSeq.R2.fq /path_to/Example.large.TellSeq.index.fq TELLSeq
```
The Minfo (optional) should be prepared as:
```
#single sample
S1=S1
S2=S2
S3=S3
S4=S4
S5=S5
#pairwise comparison
P1=S4_vs_S5
#group comparison
G1=S4,S5
G2=S1,S2,S3
```

### demo report: 
LRTK produces reproducible reports during the processing of raw sequencing data and the multiple downstream analysis.
Here is the demo report from LRTK:

<a href="https://github.com/ericcombiolab/LRTK/blob/main/demo_report/LRTK.report.html">human genome</a> 

<a href="https://github.com/ericcombiolab/LRTK/blob/main/demo_report/metagenome.multi.html">metagenome</a> 

# Troubleshooting
If any question, please feel free to contact with me or my supervisor Dr.Eric Zhang.
(Chao Yang, email: cschaoyang@comp.hkbu.edu.hk)
(Lu Zhang, email: ericluzhang@hkbu.edu.hk)
