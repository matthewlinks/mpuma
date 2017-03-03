#!/usr/bin/perl -w
use Bio::SeqIO;
use Getopt::Long;
use IPC::Open2;

my ($sffFile,$outFile,$gzip);
my @options = (
	'i=s',	\$sffFile,
	'o=s',	\$outFile,
	'z=s',	\$gzip
);
&GetOptions(@options);
die "Usage: $0 -i sff -o output.fq -z 0"
unless(defined $sffFile);
$gzip = 0 unless defined $gzip;

# setup output
my $outio;
if (defined $outFile){
	$outio = Bio::SeqIO->new(-file => ">$outFile", -format => 'fastq-illumina');
}else{
	$outio = Bio::SeqIO->new(-Fh => \*STDOUT, -format => 'fastq-illumina');
}

# open processes for reading FASTA
my $pidF = open2(\*FASTA_OUT,\*FASTA_IN,'sffinfo','-s',$sffFile) 
	or die "Error calling sffinfo -s $sffFile: $!";
close(FASTA_IN) or die "Error closing STDIN to FASTA stream: $!";
my $seqio = Bio::SeqIO->new(-fh => \*FASTA_OUT, -format => 'fasta');

# open process for reading QUAL
my $pidQ = open2(\*QUAL_OUT,\*QUAL_IN,'sffinfo','-q',$sffFile) 
	or die "Error calling sffinfo -s $sffFile: $!";
close(QUAL_IN) or die "Error closing STDIN to QUAL stream: $!";
my $qualio = Bio::SeqIO->new(-fh => \*QUAL_OUT, -format => 'qual');


while(my $seq = $seqio->next_seq()){
	my $qual = $qualio->next_seq();
	if ($seq->display_id() eq $qual->display_id()){
		my $fq = Bio::Seq::Quality->new(
			-id	=> $seq->display_id,
			-seq	=> $seq->seq,
			-qual	=> $qual->qual
		);
		$outio->write_seq($fq);
	}else{
		die "Error IDs do not line up ",$seq->display_id," vs. ",$qual->display_id;
	}
}
wait();
wait();

#undef $outio;

if (($gzip)and(defined $outFile)){
	system('gzip',$outFile);
}

exit 0;

