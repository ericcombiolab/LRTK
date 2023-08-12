#!/bin/bash

######################################################### Working Environment #################################################################
###main program
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

###output files
outFQ1="./Example.R1.fq"
outFQ2="./Example.R2.fq"
outBAM="./Example.bam"
outVCF1="./Example.small.variants.vcf"
outVCF2="./Example.large.variants.vcf"
outVCF3="./Example.small.variants.phased.vcf"

### intermediate results
#barcoded10xFQ* and nobarcoded10xFQ* could be generated from the FQCONVER command.
#The output of FQCONVER includes three pairs of FASTQ files (Example.R1/2.fq, Example.R1/2.fq.sort.wb.fq and Example.R1/2.fq.sort.wob.fq).
#The Example.R1/2.fq contains all reads.
#The Example.R1/2.fq.sort.wb.fq contains only reads with barcode annotation and is sorted based on the barcode sequence.
#The Example.R1/2.fq.sort.wob.fq contains only reads without barcode annotation and is sorted based on the genomic coordinates.

barcoded10xFQ1="./Example.R1.fq.sort.wb.fq"  
barcoded10xFQ2="./Example.R2.fq.sort.wb.fq"  
nobarcoded10xFQ1="./Example.R1.fq.sort.wob.fq"  
nobarcoded10xFQ2="./Example.R2.fq.sort.wob.fq"  

### should be prepared by users to test the WGS command
Sinfo="./sample.info"  

###database
BL10x="./database/WhiteList/white_list_10x_barcode.fa"
BLstLFR="./database/WhiteList/white_list_stlfr_barcode.fa"
BLTELLSeq="./database/WhiteList/4M-with-alts-february-2016.TELLSeq.txt"
GRCH38="./database/GRCH38/genome.fa"
SONIC="./database/sonic/GRCh38.sonic"
UNIQNESS="./database/Uniqness_map/"

######################################################### MKFQ #################################################################
###Linked-read simulator
$LRTK MKFQ -CF "./example/FQs/simulation/diploid_config" -IT 10x

######################################################### FQCONVER ##############################################################
###TENx 2 ULRF
$LRTK FQCONVER -I1 $raw10xFQ1 -I2 $raw10xFQ2 -IT 10x -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BL10x -F Yes -S Yes -T 4

######################################################## ALIGNMENT  ############################################################

$LRTK ALIGN -BQ1 $barcoded10xFQ1 -BQ2 $barcoded10xFQ2 -FQ1 $nobarcoded10xFQ1 -FQ2 $nobarcoded10xFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:Example\tSM:Example" -S Yes -M Yes -P 10x -T 4

######################################################### SNV calling ############################################################
#SNV
$LRTK SNV -B $outBAM -R $GRCH38 -A "FreeBayes" -T 4 -O $outVCF1

######################################################### SV calling #############################################################
#SV
$LRTK SV -B $outBAM -R $GRCH38 -A "Aquila" -T 4 -O $outVCF2 -V $outVCF1 -U $UNIQNESS

######################################################### Phasing ################################################################
#Phasing
$LRTK PHASE -B $outBAM -R $GRCH38 -A "HapCUT2" -T 12 -V $outVCF1 -O $outVCF3

######################################################### WGS pipeline ###########################################################
#WGS
$LRTK WGS -SI $Sinfo -OD $OUTDIR -DB $DATABASE -RG "@RG\tID:Example\tSM:Example" -T 32 

############################################################ others ##############################################################
###stLFR FQCONVER
$LRTK FQCONVER -I1 $rawstLFRFQ1 -I2 $rawstLFRFQ2 -IT stLFR -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BLstLFR -F Yes -S Yes -T 4

### stLFR alignment
barcodedstLFRFQ1="./Example.R1.fq.sort.wb.fq"
barcodedstLFRFQ2="./Example.R2.fq.sort.wb.fq"
nobarcodedstLFRFQ1="./Example.R1.fq.sort.wob.fq"
nobarcodedstLFRFQ2="./Example.R2.fq.sort.wob.fq"
$LRTK ALIGN -BQ1 $barcodedstLFRFQ1 -BQ2 $barcodedstLFRFQ2 -FQ1 $nobarcodedstLFRFQ1 -FQ2 $nobarcodedstLFRFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:Example\tSM:Example" -S Yes -M Yes -P stLFR -T 4

###TELLSeq 2 ULRF
$LRTK FQCONVER -I1 $rawTELLSeqFQ1 -I2 $rawTELLSeqFQ2 -IT TELLSeq -ID $rawTELLSeqFQi -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BLstLFR -F Yes -S Yes -T 4

barcodedtellseqFQ1="./Example.R1.fq.sort.wb.fq"
barcodedtellseqFQ2="./Example.R2.fq.sort.wb.fq"
$LRTK ALIGN -BQ1 $barcodedtellseqFQ1 -BQ2 $barcodedtellseqFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:Example\tSM:Example" -S Yes -M Yes -P TELLSeq -T 4

############################################################# end ###############################################################
