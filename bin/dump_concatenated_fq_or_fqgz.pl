#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use IPC::Open2;
my ($output);
my @options = (
	'o=s',	\$output,
);
&GetOptions(@options);

my $numFiles = @ARGV;
die "Usage: $0 -o output file2.fq file3.fastq file4.fq.gz file5.fq.gz" if ($numFiles < 1);

# redirect STDOUT to a file if it was specified
open(STDOUT,">$output") or die "Error opening $output for writting: $!" if defined $output;

# figure out how to handle the file based on file extension
foreach my $file (@ARGV){
	if(($file =~ /\.fq$/)or($file =~ /\.fastq$/)){
		 &dump_fq_compressed($file,0);
	}elsif(($file =~ /\.fq.gz$/)or($file =~ /\.fastq.gz$/)){
		&dump_fq_compressed($file,1)
	}else{
		die "Error no method to handle $file";
	}	
}
# c'est ca
exit 0;

# get Fastq. This is the fastest implementation I can think of
# this can handle gzip compression if gzip is in the path...
sub dump_fq_compressed{
	my $file = shift;
	die "Error lost file" unless defined $file;
	
	# allow for a simple signal as to the Gzip compression of the file
	my $isCompressed = shift;

	if ($isCompressed){
		open(FILE,"gzip -dc $file |") or die "Error opening $file: $!";
	}else{
		open(FILE,"<$file") or die "Error opening $file: $!";
	}

	print STDOUT while(<FILE>);
	close FILE or die "Error closing $file: $!";
}
