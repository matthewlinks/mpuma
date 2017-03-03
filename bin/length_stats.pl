#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Statistics::Descriptive;

my @lengths;
foreach my $file (@ARGV){
	my $seqio = Bio::SeqIO->new(-file => $file, -format => 'fasta');
	while(my $seq = $seqio->next_seq()){
		push(@lengths,$seq->length());
	}
}

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@lengths);
my $num = @lengths;
print join("\t",'# sequences',$num),"\n";
print join("\t",'Mean',$stat->mean()),"\n";
print join("\t",'Median',$stat->median()),"\n";
print join("\t",'Mode',$stat->mode()),"\n";
print join("\t",'Minimum',$stat->min()),"\n";
print join("\t",'Maximum',$stat->max()),"\n";
