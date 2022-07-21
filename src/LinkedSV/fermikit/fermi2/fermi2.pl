#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

# r188

die (qq/
Usage:   fermi2.pl <command> [arguments]\n
Command: unitig     generate Makefile for unitig assembly
         utglog     analyze log files generated by unitig
         mag2fmr    create FMR for multiple MAG unitig assemblies
\n/) if @ARGV == 0;

my $cmd = shift(@ARGV);
if ($cmd eq 'unitig') { &unitig(); }
elsif ($cmd eq 'mag2fmr') { &mag2fmr(); }
elsif ($cmd eq 'utglog') { &utglog(); }
else { die("ERROR: unknown command\n"); }

sub mag2fmr {
	my %opts = (l=>102, d=>3, M=>10000);
	getopts('ai:s:r:l:d:M:', \%opts);
	die (qq/fermi2.pl mag2fmr [-a] [-i in.fmr] <file1.mag.gz> [...]\n/) if @ARGV == 0;

	$opts{s} ||= gwhich("seqtk");
	$opts{r} ||= gwhich("ropebwt2");
	die ("ERROR: failed to find seqtk and ropebwt2") unless (-x $opts{s} && -x $opts{r});
	my @lines = ();
	my $prev = defined($opts{i})? $opts{i} : '';
	for my $fn (@ARGV) {
		unless (-f $fn) {
			warn("WARNING: skip non-existing file '$fn'");
			next;
		}
		my $pre = $fn =~ /(\S+)\.mag\.gz$/? $1 : $fn;
		push(@lines, qq/$pre.fmr:$fn $prev/);
		my $opt_rb2 = $prev? "-bRLi $prev" : "-bRL";
		$opt_rb2 .= " -M $opts{M}";
		my $tmp = qq/awk 'NR%2==0'|rev|tr "ACGT" "TGCA"|sort -S15G|tr "ACGT" "TGCA"|rev|$opts{r} $opt_rb2 > \$@ 2> \$@.log/;
		if (!defined($opts{a})) {
			my $genfa = qq/$opts{s} seq -nn -l0 -aq$opts{d} \$< | $opts{s} cutN -n1 - | $opts{s} seq -UL$opts{l} -l0 | gzip -1 > $pre.fa.gz/;
			my $seqs = qq/(gzip -dc $pre.fa.gz; $opts{s} seq -rl0 $pre.fa.gz)/;
			push(@lines, qq/\t$genfa; $seqs|$tmp; rm -f $pre.fa.gz/, "");
		} else {
			my $seqs = qq/($opts{s} seq -Ul0 \$<; $opts{s} seq -rUl0 \$<)/;
			push(@lines, qq/\t$seqs|$tmp/, "");
		}
		$prev = "$pre.fmr";
	}
	unshift(@lines, "all:$prev\n");

	print(join("\n", @lines), "\n");
}

sub unitig {
	my %opts = (t=>4, p=>'fmdef', l=>101, k=>-1, T=>61, o=>-1, m=>-1, s=>'100m');
	getopts('t:p:k:f:r:c:l:m:s:T:2E', \%opts);
	die (qq/
Usage:   fermi2.pl unitig [options] <in.fq>\n
Options: -p STR    output prefix [$opts{p}]
         -s STR    approximate genome size [$opts{s}]
         -2        2-pass error correction
         -l INT    primary read length [$opts{l}]
         -T INT    use INT-mer for post-trimming\/filtering [$opts{T}]
         -k INT    min overlap length during unitig construction [based on -l]
         -o INT    min overlap length during graph cleaning [based on -l]
         -m INT    min overlap length for unambiguous merging [based on -l]
         -t INT    number of threads [$opts{t}]
         -E        don't apply error correction
\n/) if (@ARGV == 0);

	die("ERROR: fermi2 doesn't work well with reads shorter than 70bp.\n") if ($opts{l} < 70);

	# get k-mer length for error correction
	my ($gsize, $k_ec1, $k_ec2);
	if ($opts{s} =~ /^([\d\.]+)([a-zA-Z])/) {
		$gsize = $1;
		if ($2 eq 'k' || $2 eq 'K') { $gsize *= 1000; }
		elsif ($2 eq 'm' || $2 eq 'M') { $gsize *= 1000000; }
		elsif ($2 eq 'g' || $2 eq 'G') { $gsize *= 1000000000; }
	} else {
		$gsize = $opts{s};
	}
	$k_ec1 = int(log($gsize) / log(2)) + 1;
	$k_ec2 = int(log($gsize) / log(2) * 1.7) + 1;
	++$k_ec1 if ($k_ec1&1) == 0;
	++$k_ec2 if ($k_ec2&1) == 0;

	$opts{k} = 51 + int(($opts{l} - 101) * .4 + .499) if $opts{k} < 0;
	$opts{m} = int($opts{l} * .75) + 1 if $opts{m} < 0;
	$opts{o} = $opts{k} + 5 if $opts{o} < 0;

	$opts{f} ||= gwhich("fermi2");
	$opts{r} ||= gwhich("ropebwt2");
	$opts{c} ||= gwhich("bfc");
	die("[E::main] failed to find the 'fermi2' executable") unless (-x $opts{f});
	die("[E::main] failed to find the 'ropebwt2' executable") unless (-x $opts{r});
	die("[E::main] failed to find the 'bfc' executable") unless (-x $opts{c});

	my @lines = ();
	push(@lines, qq/PREFIX=$opts{p}/, '');
	push(@lines, qq/EXE_FERMI2=$opts{f}/, qq/EXE_ROPEBWT2=$opts{r}/);
	push(@lines, qq/EXE_BFC=$opts{c}/, qq/GENOME_SIZE=$opts{s}/);
	push(@lines, qq/K_EC1=$k_ec1/, qq/K_EC2=$k_ec2/) if defined($opts{2});
	push(@lines, qq/K_UNITIG=$opts{k}/, qq/K_CLEAN=$opts{o}/, qq/K_TRIM=$opts{T}/, qq/K_MERGE=$opts{m}/);
	push(@lines, qq/N_THREADS=$opts{t}/, "");
	push(@lines, (-f $ARGV[0])? qq/INPUT=cat $ARGV[0]/ : qq/INPUT=$ARGV[0]/, "");
	push(@lines, "SHELL:=/bin/bash", qq/export SHELLOPTS:=errexit:pipefail/, "");

	push(@lines, qq/all:\$(PREFIX).mag.gz/, "");

	push(@lines, qq/\$(PREFIX).ec.fq.gz:/);
	if (defined $opts{E}) {
		push(@lines, (-f $ARGV[0])? qq/\tln -s $ARGV[0] \$@/ : qq/\t$(INPUT) | gzip -1 > $@/, "");
	} elsif (defined $opts{2}) {
		push(@lines, qq/\tbash -e -o pipefail -c '\$(EXE_BFC) -s \$(GENOME_SIZE)  -k \$(K_EC1) -t \$(N_THREADS) <(\$(INPUT)) <(\$(INPUT)) 2> \$@.log | gzip -1 > \$(PREFIX).ec1.fq.gz'; \\/);
		push(@lines, qq/\tbash -e -o pipefail -c '\$(EXE_BFC) -s \$(GENOME_SIZE) -Rk \$(K_EC2) -t \$(N_THREADS) <(\$(INPUT)) \$(PREFIX).ec1.fq.gz 2>> \$@.log | gzip -1 > \$\@'; \\/);
		push(@lines, qq/\trm -f \$(PREFIX).ec1.fq.gz/, "");
	} else {
		push(@lines, qq/\tbash -e -o pipefail -c '\$(EXE_BFC) -s \$(GENOME_SIZE) -t \$(N_THREADS) <(\$(INPUT)) <(\$(INPUT)) 2> \$@.log | gzip -1 > \$\@'/, "");
	}

	push(@lines, qq/\$(PREFIX).flt.fq.gz:\$(PREFIX).ec.fq.gz/);
	push(@lines, qq/\t\$(EXE_BFC) -1s \$(GENOME_SIZE) -k \$(K_TRIM) -t \$(N_THREADS) \$< 2> \$@.log | gzip -1 > \$@/, "");

	push(@lines, qq/\$(PREFIX).flt.fmd:\$(PREFIX).flt.fq.gz/);
	push(@lines, qq/\t\$(EXE_ROPEBWT2) -dNCr \$< > \$@ 2> \$@.log/, "");

	push(@lines, qq/\$(PREFIX).pre.gz:\$(PREFIX).flt.fmd/);
	push(@lines, qq/\t\$(EXE_FERMI2) assemble -l \$(K_UNITIG) -m \$(K_MERGE) -t \$(N_THREADS) \$< 2> \$@.log | gzip -1 > \$@/, "");

	push(@lines, qq/\$(PREFIX).mag.gz:\$(PREFIX).pre.gz/);
	push(@lines, qq/\t\$(EXE_FERMI2) simplify -CSo \$(K_CLEAN) -m \$(K_MERGE) -T \$(K_UNITIG) \$< 2> \$@.log | gzip -1 > \$@/, "");

	print(join("\n", @lines), "\n");
}

sub utglog {
	die("Usage: fermi2.pl utglog <prefix>\n") if @ARGV == 0;
	while (@ARGV) {
		my $pre = shift(@ARGV);
		my $fh;
		my @a = ($pre, 0, 0, 0, 0);
		open($fh, "$pre.raw.fmd.log") || die;
		while (<$fh>) {
			@a[1,2] = ($1, $2+$3+$4+$5) if /symbol counts.*\((\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)/;
		}
		close($fh);
		open($fh, "$pre.ec.fq.gz.log") || die;
		while (<$fh>) {
			$a[5] = $1 if /fmc_kmer_stat.*\s(\d+)\s+k/;
		}
		close($fh);
		open($fh, "$pre.ec.fmd.log") || die;
		while (<$fh>) {
			@a[3,4] = ($1, $2+$3+$4+$5) if /symbol counts.*\((\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)/;
		}
		close($fh);
		if (open($fh, "$pre.mag.gz.log")) {
			while (<$fh>) {
				if (/average read distance ([\d\.]+)/) {
					$a[6] = $1;
				} elsif (/approximate genome size: (\d+)/) {
					$a[7] = $1;
				}
			}
			close($fh);
		} else { $a[6] = 0; $a[7] = 0; }
		print(join("\t", @a), "\n");
	}
}

sub which
{
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file") && (-f "$x/$file");
	}
	return;
}

sub gwhich {
    my $progname = shift;
    my $addtional_path = shift if (@_);
    my $dirname = &dirname($0);
    my $tmp;

    chomp($dirname);
    if ($progname =~ /^\// && (-x $progname) && (-f $progname)) {
        return $progname;
    } elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
        return $tmp;
    } elsif (defined($dirname) && (-x "$dirname/$progname") && (-f "$dirname/$progname")) {
        return "$dirname/$progname";
    } elsif ((-x "./$progname") && (-f "./$progname")) {
        return "./$progname";
    } elsif (($tmp = &which($progname))) {
        return $tmp;
    } else {
        return;
    }
}

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
