#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($diversityFile,$fastaFile,$output,$column);
my @options = (
	'd=s',	\$diversityFile,
	'f=s',	\$fastaFile,
	'o=s',	\$output,
	'c=s',	\$column	
);
&GetOptions(@options);

usage() unless (
	defined $diversityFile
	and defined $fastaFile
	and defined $output
);
$column = 3 unless defined $column; #default is actual # reads
my $counts = parseDiversityCounts($diversityFile,$column);

my $inio = Bio::SeqIO->new(-file => "<$fastaFile", -format => 'fasta');
my $outio = Bio::SeqIO->new(-file => ">$output", -format => 'fasta');
while(my $seq = $inio->next_seq()){
	my $id = $seq->display_id();
	next unless defined $counts->{$id};
	$seq->description('Warning this file is created for input to megan and should likely NOT be used for any other use');
	for(my $i = 0; $i < $counts->{$id}; $i++){
		warn "$i for $id";
		$seq->display_id(join('_',$id,$i));
		$outio->write_seq($seq);
	}
}

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

sub usage {
	die "Usage: $0 -d diversity_file -f fasta_of_OTU_seqs -o output_fasta";
}
