#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($fasta,$wateredBLASTN,$output,$CUTOFF,$hasHeader);
my @options = (
	'i=s',	\$fasta,
	'w=s',	\$wateredBLASTN,
	'o=s',	\$output,
	'c=s',	\$CUTOFF,
	'h=s',	\$hasHeader
);
&GetOptions(@options);
usage() unless (
	defined $wateredBLASTN
	and defined $fasta
	and defined $output
);
$CUTOFF = 55 unless defined $CUTOFF;
$hasHeader = 0 unless defined $hasHeader;

my $ids2keep = parseWateredBlastn($wateredBLASTN,$CUTOFF);
my $inio = Bio::SeqIO->new(-file => "<$fasta", -format => 'fasta');
my $outio = Bio::SeqIO->new(-file => ">$output", -format => 'fasta');

while(my $seq = $inio->next_seq()){
	if (defined $ids2keep->{$seq->display_id()}){
		$outio->write_seq($seq);
	}
}

exit 0;


sub parseWateredBlastn{
	my $file = shift;
	die "Error lost wateredBLASTN file" unless defined $file;
	my $CUTOFF = shift;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	if ($hasHeader) {
		my $header = <FILE>;
		chomp $header;
	}
	my %ids2keep;
	while(<FILE>){
		chomp $_;
		my ($id,$desc,$percent,$length,$neighbour,$strand) = split(/\t/,$_);
		next if ((defined $CUTOFF)and($percent < $CUTOFF));
		next if (($percent == 0) and ($length eq 'N/A') and ($neighbour eq 'N/A'));
		$ids2keep{$id}++;
	}
	close FILE or die "Error closing $file: $!";
	return \%ids2keep;
}

sub usage{
	die "Usage: $0 -i fasta -w waterBLASTN_report -o output";
}


