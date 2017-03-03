#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my ($otuColumn,$header,$minPrevalence,$minAbundance,$abundanceCol);
my @options = (
	'o=s',	\$otuColumn,
	'h=s',	\$header,
	'p=s',	\$minPrevalence,
	'a=s',	\$minAbundance,
	'c=s',	\$abundanceCol
);
&GetOptions(@options);
$otuColumn = 0 unless defined $otuColumn;
$header = 0 unless defined $header;
$minAbundance = 0 unless defined $minAbundance;
$abundanceCol = 3 unless defined $abundanceCol; # absolute count
die "Usage: $0 -o 0 -h 0 -p 1 (N.B. prevalence is an integer for the minimum number of cases. NOT a percentage" unless defined $minPrevalence;

my %OTUprevalence;
foreach my $file (@ARGV){
	open(FILE,"<$file") or die "Error opening $file for opening: $!";
	my $head = <FILE> if ($header);
	my %seenOTUs;
	while(<FILE>){
		chomp $_;
		my @cols = split(/\s+/,$_);
		my $otu = $cols[$otuColumn];
		next if defined $seenOTUs{$otu};
		my $abundance = $cols[$abundanceCol];
		next if ($abundance < $minAbundance);
		$seenOTUs{$otu}++;
		$OTUprevalence{$otu}++;
	}
	close FILE or die "Error closing $file: $!";
}

foreach my $otu (sort {return $a cmp $b} keys %OTUprevalence){
	print $otu,"\n" if ($OTUprevalence{$otu} >= $minPrevalence);
}

exit 0;

