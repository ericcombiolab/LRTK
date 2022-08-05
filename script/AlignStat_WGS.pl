#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use File::Path;
use FindBin qw($Bin);

my ($bam, $samtools, $outdir, $move, $sam, $chr_len, $win_width);
GetOptions(
	"in:s"=>\$bam,
	"sam:s"=>\$sam,
	"samtools:s"=>\$samtools,
	"outdir:s"=>\$outdir,
	"move:s"=>\$move,
	"l:i"=>\$chr_len,
	"win_width:i" => \$win_width
);

die "Usage:perl $0 -in* <bam> -sam* <sample_name> -l <genome size> -outdir <process result dirname> -move <final result dirname> -win_width <window length,default:4000000>" unless ($bam && $sam);
$outdir ||= "./"; 
$outdir = File::Spec->rel2abs($outdir);
mkpath $outdir;
$move ||= $outdir;
$move = File::Spec->rel2abs($move);
mkpath $move;
$chr_len ||= 2897310462;
$samtools ||= "$Bin/samtools";
$win_width ||= 4000000;

my $clean_reads=0;
my $clean_bases=0;
my $mapped_reads=0;
my $mapped_bases=0;
my $mapping_rate=0;
my $dup_reads=0;
my $dup_rate=0;
my $mismatch_bases=0;
my $mismatch_rate=0;
my $uniq_reads=0;
my $uniq_bases=0;
my $uniq_rate=0;
my %insert;
open BAM,"$samtools view $bam | " or die $!;
while(<BAM>)
{
	chomp;
	my @line = split /\t/,$_;
	$insert{abs($line[8])} ++ if ( abs($line[8]) != 0 && abs($line[8]) < 700);
	my $fflag= $line[1];
        next if($fflag & 0x100);

	$clean_reads++;
	$clean_bases+=length($line[9]);
	unless ($fflag & 0x4) {   ### 0x4 presents unmap
		$mapped_reads++;
		if($_=~/XC:i:(\d+)/) {$mapped_bases+=$1;}
		else {$mapped_bases+=length($line[9]);}
		unless($fflag & 0x400)  ###0x400 presents duplicate
		{
			if($_ =~ /X0:i:(\d+)/)
			{
				$uniq_reads++ if($1==1);
				$uniq_bases+=length($line[9]);
			}
			elsif ($line[4] > 0) {
				$uniq_reads++;
				$uniq_bases+=length($line[9]);
			}
			if($_=~/XM:i:(\d+)/)
			{
				$mismatch_bases+=$1;
			}
			elsif ($_ =~ /MD:Z:(\S+)/) {
				$mismatch_bases += CoutMismathNo($_);
			}
		}
		else{
			$dup_reads++
		}
	}

}
$mismatch_rate=$mismatch_bases/$mapped_bases;
$mapping_rate=$mapped_reads/$clean_reads;
$dup_rate=$dup_reads/$mapped_reads;
$uniq_rate=$uniq_reads/$mapped_reads;

#my $name=basename($bam);
#my $sample=(split /\./,$name)[0];
open OUT, ">$move/$sam.AlignmentStat.xls" or die $!;
print OUT "Sample\t$sam\n";
print OUT "Clean reads\t$clean_reads\n";
printf OUT "Clean bases (Mb)\t%.2f\n", $clean_bases/1000000;
print OUT "Mapped reads\t$mapped_reads\n";
print OUT "Mapped bases(bp)\t$mapped_bases\n";
printf OUT "Mapping rate (%%)\t%.2f\n", 100*$mapping_rate;
print OUT "Unique reads\t$uniq_reads\n";
print OUT "Unique bases(bp)\t$uniq_bases\n";
printf OUT "Unique rate (%%)\t%.2f\n", 100*$uniq_rate;
print OUT "Duplicate reads\t$dup_reads\n";
printf OUT "Duplicate rate (%%)\t%.2f\n", 100*$dup_rate;
print OUT "Mismatch bases(bp)\t$mismatch_bases\n";
printf OUT "Mismatch rate (%%)\t%.2f\n", 100*$mismatch_rate;
close BAM;


my $basethres = 0;
my $mapQthres = 0;
my $total_chr = $chr_len;
my %depth=();
my $maxCov=0;
my $Average_sequencing_depth=0;
my $Average_sequencing_depth4=0;
my $Average_sequencing_depth10=0;
my $Average_sequencing_depth20=0;
my $Coverage=0;
my $Coverage4=0;
my $Coverage10=0;
my $Coverage20=0;
my $Coverage_bases=0;
my $Coverage_bases_4=0;
my $Coverage_bases_10=0;
my $Coverage_bases_20=0;
my $total_Coverage_bases=0;
my $total_Coverage_bases_4=0;
my $total_Coverage_bases_10=0;
my $total_Coverage_bases_20=0;
my %hashdepth;
open DEPTH, "$samtools depth  -q $basethres -Q $mapQthres $bam |" or die $!;
while(<DEPTH>){
	chomp;
	my @tmp = split /\t/,$_;
	$depth{$tmp[-1]}+=1;
	my $pos = int($tmp[1]/$win_width)+1;
	$hashdepth{$tmp[0]}{$pos}{tot} += $tmp[2];
	$hashdepth{$tmp[0]}{$pos}{num} ++;
}
close DEPTH;
my @depth=sort {$a<=>$b} keys %depth;
open HIS,">$outdir/$sam.depth_frequency.txt" or die;
open CUM,">$outdir/$sam.cumu.txt" or die;

foreach my $depth1 (sort {$a<=>$b} keys %depth){
	next if($depth1==0);
	my $per=$depth{$depth1}/$total_chr;
	$total_Coverage_bases += $depth1*$depth{$depth1};
	$Coverage_bases += $depth{$depth1};
	if($depth1>=4){
		$total_Coverage_bases_4 += $depth1 * $depth{$depth1};
		$Coverage_bases_4 += $depth{$depth1};
	}
	if($depth1>=10){
		$total_Coverage_bases_10 += $depth1 * $depth{$depth1};
		$Coverage_bases_10 += $depth{$depth1};
	}
	if($depth1>=20){
		$total_Coverage_bases_20 += $depth1 * $depth{$depth1};
		$Coverage_bases_20 += $depth{$depth1};
	}
	$maxCov=$per if($maxCov<$per);
	my $tmp=0;
	print HIS "$depth1\t$per\n";
	foreach my $depth2(@depth){
		$tmp+=$depth{$depth2} if($depth2 >= $depth1);
	}
	$tmp=$tmp/$total_chr;
	print CUM "$depth1\t$tmp\n";
}

$Average_sequencing_depth=$total_Coverage_bases/$total_chr;
$Coverage=$Coverage_bases/$total_chr;
$Average_sequencing_depth4=$total_Coverage_bases_4/$total_chr;
$Coverage4=$Coverage_bases_4/$total_chr;
$Average_sequencing_depth10=$total_Coverage_bases_10/$total_chr;
$Coverage10=$Coverage_bases_10/$total_chr;
$Average_sequencing_depth20=$total_Coverage_bases_20/$total_chr;
$Coverage20=$Coverage_bases_20/$total_chr;

print OUT "Average sequencing depth (X)\t",sprintf("%.2f",$Average_sequencing_depth),"\n";
print OUT "Coverage (\%)\t",sprintf("%.2f",100*$Coverage),"\n";
print OUT "Coverage at least 4X (\%)\t",sprintf("%.2f",100*$Coverage4),"\n";
print OUT "Coverage at least 10X (\%)\t",sprintf("%.2f",100*$Coverage10),"\n";
print OUT "Coverage at least 20X (\%)\t",sprintf("%.2f",100*$Coverage20),"\n";

close HIS;
close CUM;

#open OUT, ">$outdir/$sam.depth.window_stat.xls" or die $!;
#print OUT "chr\tpos\tdepth\n";
#foreach my $chr (sort keys %hashdepth){
#	foreach my $id (sort {$a<=>$b} keys %{$hashdepth{$chr}}){
#		$hashdepth{$chr}{$id}{num} = ($hashdepth{$chr}{$id}{num} == 0)?1:$hashdepth{$chr}{$id}{num};
#		my $mean_depth = $hashdepth{$chr}{$id}{tot} / $hashdepth{$chr}{$id}{num} ;
#		my $ids = ($id-1)*$win_width+1;
#		print OUT "$chr\t$ids\t$mean_depth\n";
#	}
#}
#close OUT;
my $ylim = 100*$maxCov;
my ($xbin,$ybin);
$ylim= int($ylim) + 1;
if($ylim <= 3){
	$ybin = 0.5;
}else{
	$ybin=1;
}
my $xlim=0;
if($Average_sequencing_depth<30){
	$xlim=100;
	$xbin=20;
}elsif($Average_sequencing_depth < 50){
	$xlim=160;
	$xbin=20;
}elsif($Average_sequencing_depth  < 120){
	$xlim=250;
	$xbin=50;
}else{
	$xlim=600;
	$xbin=100;
}
	histPlot($outdir,$move,$sam,"$outdir/$sam.depth_frequency.txt",$ylim,$ybin,$xlim,$xbin);
	cumuPlot($outdir,$move,$sam,"$outdir/$sam.cumu.txt",$xlim,$xbin);


open INSERT, ">$outdir/$sam.insert.list" or die $!;
foreach my $key( sort {$a<=>$b} keys %insert){
	print INSERT "$key\t$insert{$key}\n";
}
close INSERT;
open R, ">$outdir/$sam.Insert.plot.R" or die $!;
print R <<CODE;
rt<- read.table("$outdir/$sam.insert.list",header=F)
opar <- par()
x<-rt\$V1
y<-rt\$V2
y<-as.numeric(y)
y<-y/sum(y)*100

pdf("$move/$sam.Insert.pdf",w=8,h=6)
ymin<-0
ymax<-max(y)+0.2
ybin<-0.2
xmin<-0
xmax<-700
xbin<-100
par(mar=c(4.5, 4.5, 2.5, 2.5))
plot(x,y, ylim=c(0,ymax), xlim=c(0,xmax), col="red", type='l', lwd=1, bty="l",xaxt="n",yaxt="n", xlab="", ylab="" )

xpos <- seq(0,xmax,by=xbin)
ypos <- seq(0,ymax,by=ybin)
axis(side=1, xpos, tcl=0.2, labels=FALSE)
axis(side=2, ypos, tcl=0.2, labels=FALSE)
mtext("Insert size ($sam)", side=1, line=2, at=median(xpos), cex=1.5)
mtext("Fraction of paired reads (%)", side=2, line=3, at=median(ypos), cex=1.5 )
mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
dev.off()
CODE
close R;


my $md_value;
sub CoutMismathNo
{
	my ($align_result) = @_;
	if ($align_result =~ /MD:Z:(\S+)/){ #MD:Z:7C22C5^ACCCT46
		$md_value = $1;
	}
	else {
		print STDERR "$align_result\n";
	}

	my $mark = 0;
	my $mismatch_no = 0;
	for (my $i = 0; $i < length($md_value); $i++) {
		if (substr($md_value, $i, 1) eq '^') {
			$mark = 1;
		}
		elsif (substr($md_value, $i, 1) =~ /[A-Z]/ && $mark == 0) {
			$mismatch_no++;
		}
		elsif (substr($md_value, $i, 1) =~ /\d/) {
			$mark = 0;
		}
	}
	return $mismatch_no;
}

sub cumuPlot {
        my ($outdir, $move, $sam, $dataFile, $xlim, $xbin) = @_;
        my $figFile = "$move/$sam.Cumulative.pdf";
        my $Rline=<<Rline;
        pdf(file="$figFile",w=8,h=6)
        rt <- read.table("$dataFile",header=F)
        opar <- par()
        x <- rt\$V1[1:($xlim+1)]
        y <- 100*rt\$V2[1:($xlim+1)]
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
        par(opar)
        dev.off()
Rline
	open (ROUT,">$outdir/$sam.Cumulative.plot.R");
        print ROUT $Rline;
        close(ROUT);
}

sub histPlot {
        my ($outdir, $move, $sam, $dataFile, $ylim, $ybin, $xlim, $xbin) = @_;
        my $figFile = "$move/$sam.Depth.pdf";
        my $Rline=<<Rline;
        pdf(file="$figFile",w=8,h=6)
        rt <- read.table("$dataFile",header=F)
        opar <- par()
        t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
        y=c(rt\$V2[1:$xlim],t)
        y <- y*100
        x <- rt\$V1[1:($xlim+1)]
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
        xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Sequencing depth ($sam)",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        end <- length(xpos)-1
        mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
        par(opar)
        dev.off()
Rline
	open (ROUT,">$outdir/$sam.Depth.plot.R");
        print ROUT $Rline;
        close(ROUT);
}
