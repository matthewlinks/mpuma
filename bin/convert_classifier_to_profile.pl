#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($classifierFile,$outputFile,$term,$confidence);
my @options = (
	'c=s',	\$classifierFile,
	'o=s',	\$outputFile,
	't=s',	\$term,
	'p=s',	\$confidence
);

&GetOptions(@options);
$term = "phylum" unless defined $term;
$confidence = 0.8 unless defined $confidence;

die "Usage: $0 -c classifier -o output -t Phylum -p 0.80"
unless(
	defined $classifierFile
	and defined $term
	and defined $confidence
);

if (defined $outputFile){
	open(STDOUT,">$outputFile") or die "Error opening $outputFile for writting: $!";
}

my %results;
my $total = 0;

open(FILE,"<$classifierFile") or die "Error opening $classifierFile for reading: $!";
while(<FILE>){
	chomp $_;
	my @parts = split(/\t/,$_);
	my $otu = shift @parts;
	my $null = shift @parts;
	my ($t,$level,$p);
	while(@parts){
		$t = shift @parts;
		$level = shift @parts;
		$p = shift @parts;
		last if ($level eq $term);
	}

	my $label;
	if(($level eq $term)and($p >= $confidence)){
		$label = $t;
	}elsif($level eq $term){
		$label = 'Unknown';
	}else{
		die "Error could not find $term for $otu";
	}

	$results{$label}++;
	$total++;
}
close FILE or die "Error closing $classifierFile: $!";

print join("\t",'TaxonomicTerm','PercentComposition'),"\n";
foreach my $term (sort {return $results{$b} <=> $results{$a}} keys %results){
	my $percent = sprintf("%.2f",($results{$term}/$total)*100);
	print join("\t",$term,$percent),"\n";
}
