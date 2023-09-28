#!/bin/bash

###default
SAMTOOLS="/datahome/comp/ericteam/cschaoyang/SOFTWARE/ANACONDA/install/bin/samtools"
BEDTOOLS="/datahome/comp/ericteam/cschaoyang/SOFTWARE/ANACONDA/install/bin/bedtools"
ABD1="/datahome/comp/ericteam/cschaoyang/SOFTWARE/LRTK/LRTK2/metagenome/CalculateCoverage.round1.pl"
ABD2="/datahome/comp/ericteam/cschaoyang/SOFTWARE/LRTK/LRTK2/metagenome/CalculateCoverage.round2.pl"
ER="/datahome/comp/ericteam/cschaoyang/SOFTWARE/LRTK/LRTK2/metagenome/ExtractRef.pl"
EMA="/datahome/comp/ericteam/cschaoyang/SOFTWARE/ANACONDA/install/lib/python3.7/site-packages/script/EMA/ema"
BWA="/datahome/comp/ericteam/cschaoyang/SOFTWARE/ANACONDA/install/envs/linkedreads/bin/bwa"
RC="/datahome/comp/ericteam/cschaoyang/SOFTWARE/LRTK/LRTK2/metagenome/Readcount.Bam2bed.pl"

REF="/datahome/comp/ericteam/cschaoyang/SOFTWARE/MIDAS2/DATABASE/uhgg_merge/UHGG.reference.fa"
GENOME="/datahome/comp/ericteam/cschaoyang/SOFTWARE/MIDAS2/DATABASE/uhgg_merge/UHGG.reference.txt"
THREADS="200"

###input
#FQ1="/datahome/datasets/ericteam/cschaoyang/PROJECT/LRTK/MetaGenome/01CleanData/RESULTS/SRR14763284/SRR14763284.rmhost.1.fastq.gz"
#FQ2="/datahome/datasets/ericteam/cschaoyang/PROJECT/LRTK/MetaGenome/01CleanData/RESULTS/SRR14763284/SRR14763284.rmhost.2.fastq.gz"
#FQ1WB="/datahome/datasets/ericteam/cschaoyang/PROJECT/LRTK/MetaGenome/LRTK/MetaGenomeResults/linkedread/SRR14763284/Example.R1.fq.sort.wb.fq"
#FQ2WB="/datahome/datasets/ericteam/cschaoyang/PROJECT/LRTK/MetaGenome/LRTK/MetaGenomeResults/linkedread/SRR14763284/Example.R2.fq.sort.wb.fq"
#FQ1WOB="/datahome/datasets/ericteam/cschaoyang/PROJECT/LRTK/MetaGenome/LRTK/MetaGenomeResults/linkedread/SRR14763284/Example.R1.fq.sort.wob.fq"
#FQ2WOB="/datahome/datasets/ericteam/cschaoyang/PROJECT/LRTK/MetaGenome/LRTK/MetaGenomeResults/linkedread/SRR14763284/Example.R2.fq.sort.wob.fq"

BARCODE=$1
FQ1WB=$2
FQ2WB=$3
FQ1WOB=$4
FQ2WOB=$5
OUTDIR=$6

###output
#OUTDIR="/datahome/comp/ericteam/cschaoyang/SOFTWARE/LRTK/LRTK2/metagenome/RESULTS"
BAM1="$OUTDIR/align.round1.bam"
BAM2="$OUTDIR/align.round2.bam"

$EMA align -1 $FQ1WB -2 $FQ2WB -r $REF -R '@RG\tID:SRR14763284\tSM:SRR14763284' -p  10x  -t $THREADS | $SAMTOOLS view  -Sb  -o $OUTDIR/tmp.ema.bam --threads $THREADS
$BWA mem $REF $FQ1WOB $FQ2WOB -R '@RG\tID:SRR14763284\tSM:SRR14763284' -t $THREADS | $SAMTOOLS view  -Sb  -o $OUTDIR/tmp.bwa.bam --threads $THREADS
$SAMTOOLS merge $OUTDIR/tmp.merge.bam $OUTDIR/tmp.bwa.bam $OUTDIR/tmp.ema.bam --threads $THREADS
$SAMTOOLS sort --threads $THREADS $OUTDIR/tmp.merge.bam -o $BAM1
$SAMTOOLS index $BAM1 -@ $THREADS
$BEDTOOLS makewindows -g $GENOME -w 1000 >  $OUTDIR/windows.bed 
perl $RC -db $OUTDIR/windows.bed -bam $BAM1 -outfile $OUTDIR/windows.depth.xls
#$SAMTOOLS bedcov windows.bed $BAM > windows.depth.xls
perl $ABD1 $OUTDIR/windows.depth.xls $OUTDIR/UHGG.coverage_average.xls
rm -rf $OUTDIR/tmp.merge.bam $OUTDIR/tmp.bwa.bam $OUTDIR/tmp.ema.bam

GENOME="$OUTDIR/UHGG.reference.subset.txt"

perl $ER $OUTDIR/UHGG.coverage_average.xls $OUTDIR/UHGG.reference.subset.fa $OUTDIR/UHGG.reference.subset.txt
$BWA index $OUTDIR/UHGG.reference.subset.fa
$SAMTOOLS faidx $OUTDIR/UHGG.reference.subset.fa
$EMA align -1 $FQ1WB -2 $FQ2WB -r $OUTDIR/UHGG.reference.subset.fa -R '@RG\tID:SRR14763284\tSM:SRR14763284' -p 10x -t $THREADS | $SAMTOOLS view -Sb -o $OUTDIR/tmp.ema.bam --threads $THREADS
$BWA mem $OUTDIR/UHGG.reference.subset.fa $FQ1WOB $FQ2WOB -R '@RG\tID:SRR14763284\tSM:SRR14763284' -t $THREADS | $SAMTOOLS view -Sb -o $OUTDIR/tmp.bwa.bam --threads $THREADS
$SAMTOOLS merge $OUTDIR/tmp.merge.bam $OUTDIR/tmp.bwa.bam $OUTDIR/tmp.ema.bam --threads $THREADS
$SAMTOOLS sort --threads $THREADS $OUTDIR/tmp.merge.bam -o $BAM2
$SAMTOOLS index $BAM2 -@ $THREADS
$BEDTOOLS makewindows -g $GENOME -w 1000 > $OUTDIR/windows.iter.bed
#$SAMTOOLS bedcov $OUTDIR/windows.iter.bed $BAM2 > $OUTDIR/windows.depth.iter.xls
perl $RC -db $OUTDIR/windows.iter.bed -bam $BAM2 -outfile $OUTDIR/windows.depth.iter.xls
perl $ABD2 $OUTDIR/windows.depth.iter.xls $OUTDIR/UHGG.coverage_average.iter.xls
