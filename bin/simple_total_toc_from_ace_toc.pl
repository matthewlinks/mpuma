#!/usr/bin/perl -w
use Getopt::Long;
my ($input,$output);
my @options = (
	'i=s',	\$input,
	'o=s',	\$output
);
&GetOptions(@options);
die "Usage: $0 -i input -o output" unless (
	defined $input
	and defined $output
);


open(OUT,">$output") or die "Error opening $output for writting: $!";
print OUT "ID,nrOTU,Rationale\n";
open(IN,"<$input") or die "Error opening $input for reading: $!";
while(<IN>){
	chomp $_;
	my ($OTU,$read) = split(/\t/,$_);
	print OUT join(",",$read,$OTU,'SimpleSangerClustering'),"\n";
}
close IN or die "Error closing $input: $!";
close OUT or die "Error closing $output: $!";
exit 0;

