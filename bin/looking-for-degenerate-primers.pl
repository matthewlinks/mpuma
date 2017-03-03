#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Perl;
use Getopt::Long;

my ($nmer_length,$primer_seqs,$output,$cycleSeq);

my @options = (
	'n=s',	\$nmer_length,
	'p=s',	\$primer_seqs,
	'o=s',	\$output,
	'c=s',	\$cycleSeq
);

&GetOptions(@options);
$cycleSeq = 0 unless defined $cycleSeq;

# write to a file if specified
if (defined $output) {
	open(STDOUT,">$output") or die "Error opening $output for writting: $!";
}


$nmer_length = 26 unless defined $nmer_length;
my $INC = 100000;
#my $INC = 10000;
my $VERBOSE = 1;
my $MAXWINDOW = 30;

my $primers = read_primers($primer_seqs);

print join("\t",'FILE','SENSE Primer','RC Primer','UNKNOWN'),"\n";

foreach my $file (@ARGV){
	my %mers;

	my $sense_primers = 0;
	my $rc_primers = 0;
	my $unknown = 0;
	
	my $seqio = Bio::SeqIO->new(-file => "<$file", -format => 'fasta');
	while(my $seq = $seqio->next_seq()){
#		warn "Looking at ",$seq->display_id();
		my $match;
		my $i = 1;
		my $MIDseq = '';
		for(; $i <= $MAXWINDOW; $i++){ 
			my $end = ($nmer_length + $i) - 1;
			last if ($end > $seq->length());
			my $mer = uc($seq->subseq($i,$end));
			if (defined $primers->{$mer}){
				$match = '5primeTo3prime';
				$sense_primers++;
				my $mid_start = 1;
				my $mid_end = ($i - 1);
				if ($mid_end > 1){
					$MIDseq = $seq->subseq(1+$cycleSeq,$mid_end);
				}
			}
			last if (defined $match);
		}
		if (!defined $match){
			$i = 1;
			for(; $i <= $MAXWINDOW; $i++){
				my $end = ($seq->length() - $i)+1;
				my $start = ($end - $nmer_length)+1;
				last if ($start < 1);
				my $mer = uc($seq->subseq($start,$end));
				if (defined $primers->{$mer}){ 
	                                $match = '3primeTo5prime';
	                                $rc_primers++; 
	                                my $mid_start = $i + 1;
	                                my $mid_end = $seq->length(); 
	                                if ($mid_start < $seq->length()){ 
	                                        $MIDseq = $seq->subseq($mid_start,$mid_end);
	                                }
#					die "Mer = $mer";
	                        }elsif(defined $primers->{reverse_complement_as_string($mer)}){
					$match = '3primeTo5prime';
	                                $rc_primers++; 
	                                my $mid_start = $i + 1;
	                                my $mid_end = $seq->length(); 
	                                if ($mid_start < $seq->length()){ 
	                                        $MIDseq = reverse_complement_as_string($seq->subseq($mid_start,$mid_end));
	                                }
				}
	                	last if (defined $match); 
			}
		}

	

		$unknown++ unless defined $match;
		print join("\t",$file,$seq->display_id(),$match,$i,$MIDseq),"\n" if defined $match;
	}
#	warn("\t",$file,$sense_primers,$rc_primers,$unknown),"\n";
}

exit 0;

sub read_primers {
	my $file = shift;
	my %primers;
	my $seqio = Bio::SeqIO->new(-file => "<$file", -format => 'fasta');
	my $n = 0;
	while(my $seq = $seqio->next_seq()){
		die "Error nmer length is $nmer_length, but got a primer with a length of ",$seq->length()," ID = ",$seq->display_id()
			unless ($nmer_length == $seq->length());
		$primers{uc($seq->seq())}++;
		$n++;
		warn "Read $n" if ($VERBOSE and ($n % $INC == 0));
	}
	return \%primers;
}
