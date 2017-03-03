#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;

my ($input,$output,$blastx,$nucOutput);
my @options = (
	'i=s',	\$input,
	'o=s',	\$output,
	'b=s',	\$blastx,
	'n=s',	\$nucOutput
);
&GetOptions(@options);
usage() unless(
	defined $input
	and defined $blastx
);

my ($frames,$orientation) = parse_frames ($blastx);

my $seqio = Bio::SeqIO->new(-file => "<$input", -format => 'fasta');
my $outio;
if (defined $output){
	$outio = Bio::SeqIO->new(-file => ">$output", -format => 'fasta');
}else{
	$outio = Bio::SeqIO->new(-Fh => \*STDOUT, -format => 'fasta');
}

my $nucOutio = Bio::SeqIO->new(-file => ">$nucOutput", -format => 'fasta') if defined $nucOutput;

while(my $seq = $seqio->next_seq()){
	if (defined $frames->{$seq->display_id()}){
		# reverse complement as appropriate
		if ($orientation->{$seq->display_id()} == -1){
			$seq = $seq->revcom();
			$seq->description('Warning this sequence was reverse complemented');
		}
		$nucOutio->write_seq($seq) if defined $nucOutput; 	# write out the oriented sequence if we received a param asking for it

		# translate based on the frame we identified in the blast
		my $protein = $seq->translate(-frame => $frames->{$seq->display_id()});
		$protein->description('Warning this translation is automatically inferred from BLASTX results');
		$outio->write_seq($protein);
	}
}

exit 0;
sub usage {
	die "Usage: $0 -i input -o output -b blastx";
}

sub parse_frames {
	my $file = shift;
	die "Error lost blastx file" unless defined $file;
	my $searchio = Bio::SearchIO->new(-file => "<$file", -format => 'blast');
	my %frames;
	my %orientation;
	while ( my $result = $searchio->next_result() ) {
		my $hit = $result->next_hit;
		next unless defined $hit;
		my $hsp = $hit->next_hsp();
		next unless defined $hsp;
		die "Error saw ",$result->query_name(), " more than once" 
			if defined $frames{$result->query_name()};
			
		$frames{$result->query_name()} = $hsp->query->frame();
		$orientation{$result->query_name()} = $hsp->strand('query');


#		warn $result->query_name()," => ",$hsp->strand('query'),",",$hsp->query->frame();
	}
	return (\%frames,\%orientation);
}
