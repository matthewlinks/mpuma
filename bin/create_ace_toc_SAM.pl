#!/usr/bin/perl 
use strict;
use Getopt::Long;
my ($output,$sorted);
my @options = (
	'o=s',	\$output,
	's=s',	\$sorted,
);
&GetOptions(@options);
$sorted = 1 unless defined $sorted;

open(STDOUT,">$output") if (defined $output);

my %data;
foreach my $file (@ARGV){
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	while(<FILE>){
		chomp $_;
		next if ($_ =~ /^\@/);
		my ($read,$flag,$OTU) = split(/\t/,$_);
		next if ($OTU eq '*');
		if (!($sorted)){
			print join("\t",$OTU,$read),"\n";
		}else{
			$data{$OTU}->{$read}++;
		}
	}
	close FILE or die "Error closing $file: $!";
}

if ($sorted){
	foreach my $OTU (sort {return $a cmp $b} keys %data){
		foreach my $read (sort {return $a cmp $b} keys %{$data{$OTU}}){
			print join(' ',$OTU,$read),"\n";
		}
	}
}

exit 0;
