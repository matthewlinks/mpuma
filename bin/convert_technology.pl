#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
my ($tech,$outdir);
my @options = (
	't=s',	\$tech,
	'o=s',	\$outdir
);
&GetOptions(@options);
die "Error need technology file (-t)" unless defined $tech;
die "Error need directory to write output to (-o)" unless defined $outdir;

my @OTU_order;
open(TECH,"<$tech") or die "Error opening $tech: $!";
my $header = <TECH>;
chomp $header;
while(<TECH>){
	chomp $_;
	my ($otu) = split(/\t/,$_);
	push(@OTU_order,$otu);
}
close TECH or die "Error closing $tech: $!";


foreach my $file (@ARGV){
	
	# Read in the existing OTU information for this library
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my %OTU_data;
	while(<FILE>){
		chomp $_;
		my ($otu) = split(/\t/,$_);
		$OTU_data{$otu} = $_;
	}
	close FILE or die "Error closing $file: $!";

	my $base = basename($file);

	# write out the OTU information in EXACTLY the same order as the technology file
	my $outfile = join('/',$outdir,$base);
	open(DEST,">$outfile") or die "Error opening $outfile for writting: $!";
	print DEST $header,"\n";
	foreach my $otu (@OTU_order){
#		my $line = join("\t",$otu,'N/A','N/A','N/A','something');
		my $line = join("\t",$otu,0,0,0,'something');
		if (defined $OTU_data{$otu}){
			$line = $OTU_data{$otu};
		}
		print DEST $line,"\n";
	}
	close DEST or die "Error closing $outfile: $!";
}

exit 0;

