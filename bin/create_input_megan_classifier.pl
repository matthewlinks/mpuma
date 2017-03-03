#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($diversityFile,$classifierFile,$output,$column);
my @options = (
	'd=s',	\$diversityFile,
	'f=s',	\$classifierFile,
	'o=s',	\$output,
	'c=s',	\$column	
);
&GetOptions(@options);

# Usage
die "Usage: $0 -d diversity_file -f classifier_file -o output"
unless (
	defined $diversityFile
	and defined $classifierFile
	and defined $output
);
$column = 3 unless defined $column; #default is actual # reads
my $counts = parseDiversityCounts($diversityFile,$column);

open(OUTPUT,">$output") or die "Error opening $output for writting: $!";
open(FILE,"<$classifierFile") or die "Error opening $classifierFile: $!";
while(<FILE>){
	chomp $_;
	my @parts = split(/\t/,$_);
	my $id = shift @parts;
	next unless defined $counts->{$id};
	for(my $i = 0; $i < $counts->{$id}; $i++){
		warn "$i for $id";
		my $tmpID = join('_',$id,$i);
		print OUTPUT join("\t",$tmpID,@parts),"\n";
	}
}
close FILE;
close OUTPUT;
exit 0;

sub parseDiversityCounts {
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $col = shift;
	die "Error lost column" unless defined $col;

	my %counts;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	while(<FILE>){
		chomp $_;
		my @cols = split(/\t/,$_);
		$counts{$cols[0]} = $cols[$col];
	}
	close FILE or die "Error closing $file: $!";
	return \%counts;
}

