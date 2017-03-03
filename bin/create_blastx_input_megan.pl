#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;
use Bio::SearchIO::Writer::TextResultWriter;
use Getopt::Long;
my ($diversityFile,$output,$column,$blastx,$VERBOSE);
my @options = (
	'd=s',	\$diversityFile,
	'o=s',	\$output,
	'c=s',	\$column,
	'b=s',	\$blastx,
	'v=s',	\$VERBOSE
);
&GetOptions(@options);

usage() unless (
	defined $diversityFile
	and defined $output
	and defined $blastx
);
$column = 3 unless defined $column; #default is actual # reads
$VERBOSE = 0 unless defined $VERBOSE;


my $counts = parseDiversityCounts($diversityFile,$column);



my $inio = Bio::SearchIO->new(-file => "<$blastx", -format => 'blast');
my $writer = Bio::SearchIO::Writer::TextResultWriter->new();
my $outio = Bio::SearchIO->new(-writer => $writer , -file => ">$output");
#my $outio = Bio::SearchIO->new(-file => ">$output", -format => 'blast');
while ( my $result = $inio->next_result() ) {
	my $id = $result->query_name();
	next unless defined $counts->{$id};
	for(my $i = 0; $i < $counts->{$id}; $i++){
		warn "$i for $id" if $VERBOSE;
		$result->query_name(join('_',$id,$i));
		$outio->write_result($result);
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
	die "Usage: $0 -d diversity_file -b blastx_of_OTUs -o output -c 3";
}
