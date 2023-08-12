#!/usr/bin/perl
use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];

##cut off
my $VarQual=13;
my $DepMin=5.43;
my $DepMax=48.88;
my $AltRefRatio=0.25;

open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;
while(<IN>){
    chomp;
    if(/#/){
        print OUT $_,"\n";
    }else{
        my @a=split(/\t/,$_);
	next if($a[3] =~ /,/);
        next if($a[4] =~ /,/);
        my @b = split(";",$a[7]);
        my $cov=0;
        my ($RefDepth, $AltDepth) = (0,0);
        foreach my $ele (@b){
            my @tmp = split(/\=/,$ele);
            if($tmp[0] eq "DP"){
                $cov=$tmp[1];
            }
            if($tmp[0] eq "RO"){
                $RefDepth = $tmp[1];
            }
            if($tmp[0] eq "AO"){
                $AltDepth = (split(",",$tmp[1]))[0];
            }
        }
        next if($cov < $DepMin);
        next if($AltDepth < 2);
	

        my $ratio=$AltDepth/($RefDepth+$AltDepth);
        if($ratio > 0.15){
            print OUT $_,"\n";
        }
    }
}
close IN;
close OUT;

