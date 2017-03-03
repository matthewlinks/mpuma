#!/usr/bin/perl -w
use strict;
use Getopt::Long;

$| = 1;

my ($input_file,$output_file,$header);
my @options = (
	'i=s',	\$input_file,
	'o=s',	\$output_file,
	'h=s',	\$header,
);
&GetOptions(@options);
usage() unless (
	defined $input_file
);
$header = 1 unless defined $header;

# redirect stdout if there is a file to write to
if (defined $output_file){
	open(STDOUT,">$output_file") or die "Error opening $output_file for writting: $!" if (defined $output_file);
}

open(FILE,"<$input_file") or die "Error opening $input_file for reading: $!";
my $header_line = <FILE> if ($header);
my %possibleChimeras;
while(<FILE>){
	chomp $_;
	my @parts = split(/\t/,$_);
	# skip the IDs with no obvious concerns
	next if (($parts[1] =~ /^TRUE$/i)and($parts[2] =~ /^TRUE$/i)); 
	$possibleChimeras{$parts[0]}++;
}
close FILE;

foreach my $ID (sort {return $a cmp $b} keys %possibleChimeras){
	print $ID,"\n";
}

exit 0;
