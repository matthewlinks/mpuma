#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($config,$midString,$file,$midFile,$only_singles);

my @options = (
	'c=s',	\$config,
	'm=s',	\$midString,
	'mf=s',	\$midFile,
	'i=s',	\$file,
	'o=s',	\$only_singles,
);
&GetOptions(@options);

usage() unless defined $file
	and ((defined $midString)xor(defined $midFile));

$config = '/usr/local/rig/config/MIDConfig.parse' unless defined $config;
my $midDetails = read_config($config);
$only_singles = 1 unless defined $only_singles; # this 


# get the list of mids to parse for 
my @mids;
if (defined $midString){
	@mids = split(/\,/,$midString);
}elsif(defined$midFile){
	open(MIDLIST,"<$midFile") or die "Error opening $midFile for reading: $!";
	while(<MIDLIST>){
		chomp $_;
		my ($m) = split(/\s+/,$_);
		push(@mids,$m);
	}
	close MIDLIST or die "Error closing $midFile: $!";
}else{
	die "Error no mid file defined (-mf) and no mid string (-m)";
}

print join("\t",'Number of Matches','Read','Direction of UT primer','distance from key to primer','Matching mids'),"\n";

open(FILE,"<$file") or die "Error opening $file: $!";
while(<FILE>){
	chomp $_;
	my ($file,$read,$direction,$distance,$extra) = split(/\t/,$_);
	my $possible_mid_seq;
	if ($extra =~ /^([AGCT]+)/){
		$possible_mid_seq = $1;
	}
	next unless defined $possible_mid_seq;
	my $matching_mids = search_mids(\@mids,$midDetails,$possible_mid_seq);
	my $num_matches = @{$matching_mids};
	if ($only_singles){
		print join("\t",$num_matches,$read,$direction,$distance,$possible_mid_seq,join(',',@{$matching_mids})),"\n" if ($num_matches == 1);
	}else{
		print join("\t",$num_matches,$read,$direction,$distance,$possible_mid_seq,join(',',@{$matching_mids})),"\n";
	}
}

exit 0;

sub usage {
	die "Usage: $0 -c configOfMids -i file -m mid1,mid2 -mf midFile";
}

sub search_mids {
	my $mids = shift;
	die "Error lost mids" unless defined $mids;
	my $midDetails = shift;
	die "Error lost midDetails" unless defined $midDetails;
	my $needle = shift;
	die "Error lost needle" unless defined $needle;
	my $length = length($needle);
	
	my @matches;
	foreach my $mid (@{$mids}){
		my $seq = $midDetails->{$mid}->{seq};
		if ($seq =~ /$needle$/){
			push(@matches,$mid);
		}
	}
	return \@matches;
}

sub read_config {
	my $config = shift;
	die "Error lost file to parse" unless defined $config;
	open(FILE,"<$config") or die "Error opening $config: $!";
	my %midDetails;
	while(<FILE>){
		chomp $_;
		if ($_ =~ /\s*mid\s+\=\s+\"(\w+)\",\s*\"(\w+)\",\s*(\d+);/){
			my $name = $1;
			$midDetails{$name}->{seq} = $2;	
			$midDetails{$name}->{errors} = $3;	
		}
	}
	close FILE or die "Error closing $config: $!";
	return \%midDetails;
}
