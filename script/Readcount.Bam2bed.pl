#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($dbfile, $inbam, $outfile, $help);
GetOptions(
        "db=s"                 => \$dbfile,
        "bam=s"                => \$inbam,
        "outfile=s"            => \$outfile,
	"h|?"                  => \$help,
);

my $usage = <<USAGE;
Description: To calcualte the gene abundance based on IGC
Usage: perl $0  [options] 
Options:
	-db         <path>      IGC annotation file
	-bam        <path>      bam file
        -outfile    <path>      output file
	-h|?:  print help information

e.g.: perl $0 -infile  -outdir example/
USAGE

die $usage if (!defined $dbfile || !defined $inbam || !defined $outfile || $help);
##Gene abundance units: single read

my (%info, %uniq, %MultipleReadsAbd, %MultipleGeneAbd, %MultipleSign);
my ($AllInsertSize, $AverageInsertSize, $PairNumber) = (0, 0, 0);

###process bed
open DB,"<$dbfile" or die $!;
while(<DB>){
        chomp;
        next if(/^ID/);
        my @a = split(/\s+/,$_);
        $info{$a[0]}{$a[1]} = $a[2] - $a[1];
}
close DB;

###Processing Bamfile
#http://picard.sourceforge.net/explain-flags.html
open BAM,"samtools view $inbam |" or die $!;
while(<BAM>){
	chomp;
	next if(/^@/);
	my @a = split(/\t/,$_);
	my @flag = split(//,sprintf("%b", $a[1]));
	my $count = scalar(@flag);
	for(my $i=0; $i < 12 - $count; $i++){
		unshift @flag, 0;
	}
	#print @flag,"\n";
	my $readid= $a[0]."".$flag[-7];
	if($flag[-9] == 1 or $flag[-12] == 1){
		if(exists $MultipleSign{$readid}){
			#		print $readid,"\n";
			$MultipleSign{$readid}++;
		}else{
			$MultipleSign{$readid} = 1;
			$MultipleSign{$readid}++;
		}
	}
}
close BAM;

open BAM,"samtools view $inbam |" or die $!;
while(<BAM>){
        chomp;
        next if(/^@/);
        my @a = split(/\t/,$_);
        my $genome = $a[2];
	next if($genome eq "*");
	my $pos = $a[3];
	
	my $tmp = $info{$genome};
	my $pos_key = 0;
	for my $ele (keys %$tmp){
		if($ele < $pos and ($pos - $ele) <= $info{$genome}{$ele}){
			$pos_key = $ele;
		}
	}
	#print join("\t","###",$genome, $pos_key, $info{$genome}{$pos_key}),"\n";
	my $gname = $genome."\t".$pos_key."\t".$info{$genome}{$pos_key};

        my @flag = split(//,sprintf("%b", $a[1]));
        my $count = scalar(@flag);
        for(my $i=0; $i < 12 - $count; $i++){
                unshift @flag, 0;
        }

	next if($flag[-3]==1);
	#next if($flag[-3]==1 or $flag[-12]==1); ###To remove unmapped or supplementary aligned segment
        my $readid= $a[0]."".$flag[-7];
	if(!exists $MultipleSign{$readid}){ ###unique mapping
		$uniq{$gname} = (exists $uniq{$gname})? ($uniq{$gname}) + 1 : 1;
	}else{
                $MultipleReadsAbd{$readid} = (exists $MultipleReadsAbd{$readid})? $MultipleReadsAbd{$readid}.",".$gname : $gname;
        }
}

###Calculate the Co value
foreach my $readid (keys %MultipleReadsAbd){
	#print $readid,"\n";
        my @a = split(/\,/,$MultipleReadsAbd{$readid}); 
        my $sum = 0;
        foreach my $gname (@a){
                if(exists $uniq{$gname}){
                        $sum += $uniq{$gname};  
                }
        }
        next if($sum == 0);
        foreach my $gname (@a){
                next if(! exists $uniq{$gname});
                $MultipleGeneAbd{$gname} = (exists $MultipleGeneAbd{$gname})? $MultipleGeneAbd{$gname}+$uniq{$gname}/$sum : $uniq{$gname}/$sum;
        }
}

open OUT,">$outfile" or die $!;
print OUT join("\t","#GeneID","Start","Windows","GeneLength","UniqueReads","MultipleReads","TotalReads","Relative_abundance"),"\n";
my $sum = 0;
foreach my $gname (sort keys %uniq){
        $MultipleGeneAbd{$gname} = 0 if(! exists $MultipleGeneAbd{$gname});
        my $total = $uniq{$gname} + $MultipleGeneAbd{$gname};
	my @tmp = split(/\t/,$gname);
        my $abd   = $total/$info{$tmp[0]}{$tmp[1]};
        $sum += $abd;
}

foreach my $gname (sort keys %uniq){
        next if($uniq{$gname} == 0);
	my @tmp = split(/\t/,$gname);
        my $glen  = $info{$tmp[0]}{$tmp[1]};
        my $total = $uniq{$gname} + $MultipleGeneAbd{$gname};
        my $RelativeAbd   = $total/($glen*$sum);
        print OUT join("\t",$gname,$glen,$uniq{$gname},$MultipleGeneAbd{$gname},$total,$RelativeAbd),"\n";
}
close OUT;
