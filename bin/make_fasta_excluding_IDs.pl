#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

my ($input_file,$output_file,$exclude_file);
my @options = (
	'i=s',	\$input_file,
	'o=s',	\$output_file,
	'e=s',	\$exclude_file,
);
&GetOptions(@options);
die "Usage: $0 -i input -e excludIDfile" unless (
	defined $input_file
);

my $excludeIDs = get_IDs($exclude_file);

my $outio;
if (defined $output_file){
	$outio = Bio::SeqIO->new(-file =>">$output_file", -format => 'fasta');
}else{
	$outio = Bio::SeqIO->new(-Fh => \*STDOUT, -format => 'fasta');
}

my $inio = Bio::SeqIO->new(-file => "<$input_file", -format => 'fasta');
while(my $seq = $inio->next_seq()){
	next if (defined $excludeIDs->{$seq->display_id()});
	$outio->write_seq($seq);
}

exit 0;

sub get_IDs{
	my $file = shift;
	die "Error lost file" unless defined $file;
	my %IDs;
	open(FILE,"<$file") or die "Error opening $file for reading : $!";
	while(<FILE>){
		chomp $_;
		my ($id) = split(/\s+/,$_);
		$IDs{$id}++;
	}
	close FILE or die "Error closing $file: $!";
	return \%IDs;
}
	
