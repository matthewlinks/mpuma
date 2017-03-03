#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

my ($inputFasta,$fivePrime,$threePrime,$header,$output,$wateredBlast);
my @options = (
	'i=s',	\$inputFasta,
	'h=s',	\$header,
	'f=s',	\$fivePrime,
	't=s',	\$threePrime,
	'o=s',	\$output,
	'w=s',	\$wateredBlast
);
&GetOptions(@options);

# defaults
$fivePrime = '-5prime' unless defined $fivePrime;
$threePrime = '-3prime' unless defined $threePrime;
$header = 1 unless defined $header;

die "Usage: $0 -i inputFasta -h headerBoolean -f \'-5prime\' -t \'-3prime' -o output"
	unless (
		defined $inputFasta
		and defined $wateredBlast
		and defined $fivePrime
		and defined $threePrime
	);

# Create an odered list of the sequences
my @orderedSequences;
my $seqio = Bio::SeqIO->new(-file => "<$inputFasta", -format => "fasta");
while(my $seq = $seqio->next_seq()){
	push(@orderedSequences,$seq->display_id());
}

# parse the wateredBlast data
my %results;
open(FILE,"<$wateredBlast") or die "Error opening $wateredBlast for reading: $!";
if ($header){
	my $line = <FILE>;
}
while(<FILE>){
	chomp $_;
	my @parts = split(/\t/,$_);
	if ($parts[0] =~ /^(.*)$fivePrime$/){
		my $id = $1;
		# line is for the 5' end of the sequence
		die "Error multiples defined for ",$parts[0] if defined $results{$id}->{'5prime'};
		$results{$id}->{'5prime'}->{hit} = $parts[4];
		$results{$id}->{'5prime'}->{strand} = $parts[5];
		$results{$id}->{'5prime'}->{identity} = $parts[2];

	}elsif($parts[0] =~ /^(.*)$threePrime$/){
		my $id = $1;
		# line is for the 5' end of the sequence
		die "Error multiples defined for $parts[0]" if defined $results{$id}->{'3prime'};
		$results{$id}->{'3prime'}->{hit} = $parts[4];
		$results{$id}->{'3prime'}->{strand} = $parts[5];
		$results{$id}->{'3prime'}->{identity} = $parts[2];	
	}else{
		die "Error got unparsable ID in ",$parts[0];
	}
}
close FILE or die "Error closing $wateredBlast: $!";

# print the results
open(STDOUT,">$output") if defined $output;
print join("\t",'ID','Match Test','Strand Test','% Identity difference','5prime hit','5prime % identity','5prime strand','3prime hit','3prime % identity', '3prime strand'),"\n";
foreach my $id (@orderedSequences){
	my ($strandTest,$matchTest) = ('N/A','N/A');
	if ((defined $results{$id}->{'3prime'})and(defined $results{$id}->{'5prime'})){
		$strandTest = 'FALSE';
		$strandTest = 'TRUE' 
			if (
				(
					($results{$id}->{'3prime'}->{strand} eq '+')
					and
					($results{$id}->{'5prime'}->{strand} eq '-')
				)
				or
				(
					($results{$id}->{'3prime'}->{strand} eq '-')
					and
					($results{$id}->{'5prime'}->{strand} eq '+')
				)
			);
		$matchTest = 'FALSE';
		$matchTest = 'TRUE' if ($results{$id}->{'5prime'}->{hit} eq $results{$id}->{'3prime'}->{hit});
	}
	my $percentIdentityDifference = 'N/A';
	if ((defined $results{$id}->{'3prime'}->{identity})and(defined $results{$id}->{'5prime'}->{identity})){
		$percentIdentityDifference = sprintf("%.2f",abs($results{$id}->{'3prime'}->{identity} - $results{$id}->{'5prime'}->{identity}));
	}
	print join("\t",$id,$matchTest,$strandTest,$percentIdentityDifference,
		$results{$id}->{'5prime'}->{hit},$results{$id}->{'5prime'}->{identity},$results{$id}->{'5prime'}->{strand},
		$results{$id}->{'3prime'}->{hit},$results{$id}->{'3prime'}->{identity},$results{$id}->{'3prime'}->{strand}),"\n";
}

exit 0;


