#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use IO::File;

use Getopt::Long;
my ($VERBOSE,$fasta_file,$n_chunks,$prefix,$outdir);
my @options = (
        'i=s',  \$fasta_file,
        'c=s',  \$n_chunks,
        'p=s',  \$prefix,
        'o=s',  \$outdir,
        'v=s',  \$VERBOSE
);
&GetOptions(@options);
usage() unless defined $fasta_file;
usage() unless defined $n_chunks;
usage() unless defined $prefix;
usage() unless defined $outdir;

split_into_chunks($fasta_file,$n_chunks,$prefix,$outdir);
exit 0;

sub split_into_chunks {
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $n_chunks = shift;
	die "Error lost the number of chunks" unless defined $n_chunks;
	my $prefix = shift;
	die "Error lost prefix" unless defined $prefix;
	my $outdir = shift;
	die "Error lost outdir" unless defined $outdir;

	my $outputs = get_files($outdir,$prefix,$n_chunks);
	die "Error lost output streams" unless defined $outputs;

	my $inio = Bio::SeqIO->new(-file => "<$file", -format => 'largefasta');
	
	my $i = 0;
	while(my $seq = $inio->next_seq()){
		$outputs->[$i]->write_seq($seq);
		$i++;
		$i = 0 unless ($i < @{$outputs});
	}
}

sub get_files {
	my $dir = shift;
	die "Error lost dir" unless defined $dir;
	my $prefix = shift;
	die "Error lost prefix" unless defined $prefix;
	my $num = shift;
	die "Error lost num" unless defined $num;

	my @FHs;

	for(my $i = 1; $i <= $num; $i++){
		my $file = $dir . '/' . $prefix . '.' . $i;
		die "Error $file exists" if (-e $file);
#		my $seqio = Bio::SeqIO->new(-file => ">$file", -format => 'largefasta');
		my $seqio = Bio::SeqIO->new(-file => ">$file", -format => 'fasta');
		push(@FHs,$seqio);
	}
	return \@FHs;
}


sub usage {
	die <<EOF
Usage: $0 -i fasta_to_split -c num_chunks -p output_file_prefix -o outdir
Options:
	-v 1 verbose
EOF
}

