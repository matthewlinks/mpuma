#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;
use Bio::SearchIO::Writer::TextResultWriter;
use Getopt::Long;
use IO::Zlib;
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

my %bitscores;
while (my $result = $inio->next_result() ) {
	my $id = $result->query_name();
	warn "Reading BLAST data for $id" if $VERBOSE;
	if (defined $counts->{$id}){
		# then we need to parse these hits
		my $num_hits = 0;
		while(my $hit = $result->next_hit()){
			# hack 29Nov2013 to loose the ID which seemed to through off megan
			my ($acc,@parts) = split(/\s+/,$hit->description());
			my $meganTxt;
			if ($parts[0] =~ /^[a-z]/) {
				# then this can't start with Genus...
				$meganTxt = join(' ',$acc,@parts);
			}else{
				$meganTxt = join(' ',@parts);
			}
			$meganTxt =~ s/\,/ /g; # remove commas
#			$bitscores{$id}->{$hit->description()} = $hit->bits();
			$bitscores{$id}->{$meganTxt} = $hit->bits();
			$num_hits++;
		}
		$bitscores{$id}->{'Unknown'} = 0 if ($num_hits == 0);
	} #else skip those hits...
}

# Open the output file as a Gzip stream
my $fh = IO::Zlib->new("$output","wb9");
die "Error opening gzip file $output: $!" unless defined $fh;

foreach my $id (keys %bitscores){
	next unless defined $counts->{$id};
	for(my $i = 0; $i < $counts->{$id};$i++){
		foreach my $hit (keys %{$bitscores{$id}}){
			my $fake_read = $id . '_' . $i;
			print $fh join(',',$id . '_' . $i,$hit,$bitscores{$id}->{$hit}),"\n";
		}
	}
}
$fh->close or die "Error closing $output: $!";
undef $fh;

exit 0;

sub parseDiversityCounts {
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $col = shift;
	die "Error lost column" unless defined $col;

	my %counts;
	warn "Reading diversity for $file" if $VERBOSE;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	while(<FILE>){
		chomp $_;
		my @cols = split(/\t/,$_);
		$counts{$cols[0]} = $cols[$col];
		warn "Saw count for ",$cols[0]," => ",$cols[$col] if $VERBOSE;
	}
	close FILE or die "Error closing $file: $!";
	return \%counts;
}

sub usage {
	die "Usage: $0 -d diversity_file -b blastx_of_OTUs -o output.gz -c 3";
}
