#!/bin/bash
###INPUT && OUTPUT
export C_INCLUDE_PATH="/home/comp/zmzhang/software/htslib_install/include":$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH="/home/comp/zmzhang/software/htslib_install/include":$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH="/home/comp/zmzhang/software/htslib_install/lib":$LD_LIBRARY_PATH
export LIBRARY_PATH="/home/comp/zmzhang/software/htslib_install/lib":$LIBRARY_PATH
export PATH="/home/comp/zmzhang/software/htslib_install/bin":$PATH

curP=`pwd`
curP="/tmp/local/cschaoyang/SOFTWARE/LRTK/example"
curP="./"
raw10xFQ1=${curP}"/FQs/Example.10x.R1.fq"
raw10xFQ2=${curP}"/FQs/Example.10x.R2.fq"
rawstLFRFQ1=${curP}"/FQs/Example.stLFR.R1.fq"
rawstLFRFQ2=${curP}"/FQs/Example.stLFR.R2.fq"
rawTELLSeqFQ1=${curP}"/FQs/Example.TELLSeq.R1.fq"
rawTELLSeqFQ2=${curP}"/FQs/Example.TELLSeq.R2.fq"
rawTELLSeqFQi=${curP}"/FQs/Example.TELLSeq.index.fq"

###OUTFILE
outFQ1=${curP}"/test/Example.R1.fq"
outFQ2=${curP}"/test/Example.R2.fq"
outBAM=${curP}"/test/Example.bam"
outVCF1=${curP}"/test/Example.small.variants.vcf"
outVCF2=${curP}"/test/Example.large.variants.vcf"
outVCF3=${curP}"/test/Example.small.variants.phased.vcf"

###DATABASE
DATABASE=${curP}"/../database"
BL10x=${DATABASE}"/WhiteList/white_list_10x_barcode.fa"
BLstLFR=${DATABASE}"/WhiteList/white_list_stlfr_barcode.fa"
BLTELLSeq=${DATABASE}"/WhiteList/4M-with-alts-february-2016.TELLSeq.txt"
GRCH38=${DATABASE}"/GRCH38/genome.fa"
SONIC=${DATABASE}"/sonic/GRCh38.sonic"
UNIQNESS=${DATABASE}"/Uniqness_map/"

###main program
LRTK=${curP}"/../LRTK.py"

######################################################### MKFQ #################################################################
#Linked-read simulator
#python $LRTK MKFQ -CF ${curP}"/FQs/simulation/diploid_config" -IT 10x

######################################################### FQCONVER ##############################################################
###TENx 2 ULRF
#python $LRTK FQCONVER -I1 $raw10xFQ1 -I2 $raw10xFQ2 -IT 10x -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BL10x -F Yes -S Yes -T 4

######################################################## ALIGNMENT  ############################################################
barcoded10xFQ1=${curP}"/test/Example.R1.fq.sort.wb.fq"
barcoded10xFQ2=${curP}"/test/Example.R2.fq.sort.wb.fq"
nobarcoded10xFQ1=${curP}"/test/Example.R1.fq.sort.wob.fq"
nobarcoded10xFQ2=${curP}"/test/Example.R2.fq.sort.wob.fq"
#python $LRTK ALIGN -BQ1 $barcoded10xFQ1 -BQ2 $barcoded10xFQ2 -FQ1 $nobarcoded10xFQ1 -FQ2 $nobarcoded10xFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:NA12878\tSM:NA12878" -S Yes -M Yes -P 10x -T 4

######################################################### SNV calling ############################################################
#SNV
#python $LRTK SNV -B $outBAM -R $GRCH38 -A "FreeBayes" -T 4 -O $outVCF1

######################################################### SV calling #############################################################
#SV
#python $LRTK SV -B $outBAM -R $GRCH38 -A "Aquila" -T 4 -O $outVCF_Aquila -V $outVCF1 -U $UNIQNESS

######################################################### Phasing ################################################################
#Phasing
#python $LRTK PHASE -B $outBAM -R $GRCH38 -A "HapCUT2" -T 12 -V $outVCF1 -O $outVCF3

######################################################### WGS pipeline ###########################################################
#WGS
si=${curP}"/sample.info"
od=${curP}"/WGS/"
python $LRTK WGS -SI $si -OD $od -DB $DATABASE -RG "@RG\tID:NA12878\tSM:NA12878" -T 32

############################################################ others ##############################################################
###stLFR FQCONVER
#python $LRTK FQCONVER -I1 $rawstLFRFQ1 -I2 $rawstLFRFQ2 -IT stLFR -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BLstLFR -F Yes -S Yes -T 4

### stLFR alignment
#barcodedstLFRFQ1=${curP}"/test/Example.R1.fq.sort.wb.fq"
#barcodedstLFRFQ2=${curP}"/test/Example.R2.fq.sort.wb.fq"
#nobarcodedstLFRFQ1=${curP}"/test/Example.R1.fq.sort.wob.fq"
#nobarcodedstLFRFQ2=${curP}"/test/Example.R2.fq.sort.wob.fq"
python $LRTK ALIGN -BQ1 $barcodedstLFRFQ1 -BQ2 $barcodedstLFRFQ2 -FQ1 $nobarcodedstLFRFQ1 -FQ2 $nobarcodedstLFRFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:foo\tSM:bar" -S Yes -M Yes -P stLFR -T 4

###TELLSeq 2 ULRF
#python $LRTK FQCONVER -I1 $rawTELLSeqFQ1 -I2 $rawTELLSeqFQ2 -IT TELLSeq -ID $rawTELLSeqFQi -O1 $outFQ1 -O2 $outFQ2 -OT ULRF -B $BLstLFR -F Yes -S Yes -T 4

#barcodedtellseqFQ1=${curP}"/test/Example.R1.fq.sort.wb.fq"
#barcodedtellseqFQ2=${curP}"/test/Example.R2.fq.sort.wb.fq"
python $LRTK ALIGN -BQ1 $barcodedtellseqFQ1 -BQ2 $barcodedtellseqFQ2 -R $GRCH38 -O $outBAM -RG "@RG\tID:fo\tSM:bar" -S Yes -M Yes -P TELLSeq -T 4

############################################################# end ###############################################################
