#!/usr/bin/perl -w
use strict;
use Bio::Perl;
use Bio::SeqIO;
use Getopt::Long;
my ($fasta,$wateredBLASTN,$output);
my @options = (
	'i=s',	\$fasta,
	'w=s',	\$wateredBLASTN,
	'o=s',	\$output);
&GetOptions(@options);
usage() unless (
	defined $wateredBLASTN
	and defined $fasta
	and defined $output
);


my $strandedness = parseWateredBlastnStrandedness($wateredBLASTN);

my $inio = Bio::SeqIO->new(-file => "<$fasta", -format => 'fasta');
my $outio = Bio::SeqIO->new(-file => ">$output", -format => 'fasta');

while(my $seq = $inio->next_seq()){
	if (defined $strandedness->{$seq->display_id()}){
		my @strands = keys %{$strandedness->{$seq->display_id()}};
		my $num = @strands;
		if ($num == 1){
			my $dominant = shift @strands;
			if ($dominant eq '+'){
				$seq->description('Sequence has not been altered');
			}elsif($dominant eq '-'){
				$seq->description('Warning this sequence has been reverse complemented');
				$seq->seq(reverse_complement_as_string($seq->seq()));
			}else{
				die "Error got a strand of [$dominant] for ",$seq->display_id();
			}
		}elsif($num > 1){
			# multiple different directions detected
			warn $seq->display_id(),": has multiple strand matches so exlcuding it from the output";
			next;
		}else{
			die "Error lost strand info for ",$seq->display_id();
		}
	}else{
		$seq->description('No hit found');
	}
	$outio->write_seq($seq);
}

exit 0;


sub parseWateredBlastnStrandedness{
	my $file = shift;
	die "Error lost wateredBLASTN file" unless defined $file;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $header = <FILE>;
	chomp $header;
	my %strandedness;
	while(<FILE>){
		chomp $_;
		my ($id,$desc,$percent,$length,$neighbour,$strand) = split(/\t/,$_);
		next if (($percent == 0) and ($length eq 'N/A') and ($neighbour eq 'N/A'));
		$strandedness{$id}->{$strand}++;
	}
	close FILE or die "Error closing $file: $!";
	return \%strandedness;
}

sub usage{
	die "Usage: $0 -i fasta -w waterBLASTN_report -o output";
}


