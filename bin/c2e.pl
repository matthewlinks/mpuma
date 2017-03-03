#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
my ($inputFasta,$outputFasta,$length,$percent);

my @options = (
	'p=s',	\$percent,
	'l=s',	\$length,
	'o=s',	\$outputFasta,
	'i=s',	\$inputFasta
);
&GetOptions(@options);

die "Usage: $0 -i inputFasta -o outputFasta -p 25 (xor -l 150)"
	unless (
		defined $inputFasta
		and defined $outputFasta
		and ((defined $percent)xor(defined $length))
	);

my $seqio = Bio::SeqIO->new(-file => "<$inputFasta", -format => "fasta");
my $outio = Bio::SeqIO->new(-file => ">$outputFasta", -format => "fasta");

while(my $seq = $seqio->next_seq()){
	my ($fivePrimeSeq,$threePrimeSeq);
	my ($start,$stop) = (1,$seq->length());
	if (defined $length){
		if ($length < $seq->length){
			$stop = $length;
		}
	}elsif(defined $percent){
		my ($start,$stop) = (1,$seq->length());
		if ($percent < 100){
			$stop = int(($percent / 100)*$seq->length());
		}
	}
	$fivePrimeSeq = Bio::PrimarySeq->new(
		-seq => $seq->subseq($start,$stop),
		-display_id	=> $seq->display_id() . '-5prime'
	);
	$outio->write_seq($fivePrimeSeq);
	$threePrimeSeq = Bio::PrimarySeq->new(
		-seq => $seq->revcom->subseq($start,$stop),
		-display_id	=> $seq->display_id() . '-3prime'
	);
	$outio->write_seq($threePrimeSeq);
}



