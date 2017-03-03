#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

my ($input,$output,$mask_file);
my @options = (
	'i=s',	\$input,
	'o=s',	\$output,
	'm=s',	\$mask_file
);
&GetOptions(@options);

usage() unless (
	defined $mask_file
);

my $pos = parse_masking($mask_file);
my $inio;
if (defined $input){
	$inio = Bio::SeqIO->new(-file => "<$input", -format => 'fasta');
}else{
	$inio = Bio::SeqIO->new(-Fh => \*STDIN, -format => 'fasta');
}

my $outio;
if (defined $output){
	$outio = Bio::SeqIO->new(-file => ">$output", -format => 'fasta');
}else{
	$outio = Bio::SeqIO->new(-Fh => \*STDOUT, -format => 'fasta');
}

while(my $seq = $inio->next_seq()){
	warn "Masking ",$seq->display_id();
	foreach my $start (keys %{$pos->{$seq->display_id}}){
		die "Error $start on ",$seq->display_id," is not numeric" unless ($start =~ /^\d+$/);
		foreach my $end (keys %{$pos->{$seq->display_id}->{$start}}){
			die "Error $end on ",$seq->display_id," is not numeric" unless ($end =~ /^\d+$/);

			die "Error start is not <= end: $start $end" unless ($start <= $end);
			my $length_mask = ($end - $start)+1;
			my $ns = 'N' x $length_mask;
			
			my $masked_sequence;
			
			# if we mask after position 1 then we have to add the prefix
			if ($start > 1){
				$masked_sequence .= $seq->subseq(1,($start-1));
			}
	
			$masked_sequence .=  'N' x $length_mask;
	
			# if we stop masking at least 1bp before the end of the sequence we have to add the suffix
			if ($end < $seq->length()){
				$masked_sequence .= $seq->subseq(($end +1),$seq->length());
			}
			$seq->seq($masked_sequence);
		}
	}
	$outio->write_seq($seq);
}
	

exit;



sub parse_masking {
	my $file = shift;
	die "Error lost masking file" unless defined $file;
	my %pos;
	open(FILE,"<$file") or die "Error opening $file: $!";
	while(<FILE>){
		chomp $_;
		my ($id,$start,$end) = split(/\s+/,$_);
		$pos{$id}->{$start}->{$end}++;
	}
	close FILE or die "Error closing $file: $!";
	return \%pos;
}
sub usage {
	die "Usage: $0 -i input -o output -m filename";
}
