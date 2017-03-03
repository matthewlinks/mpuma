#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use GD::Graph::histogram;
my ($file, $output,$binsize);
my @options = (
	'i=s',	\$file,
	'o=s',	\$output,
	'b=s',	\$binsize,
);
GetOptions(@options);

$binsize = 5 unless defined $binsize;


# get the max % identity for each record
my %best_ids;
open(FILE,"<$file") or die "Error opening $file for reading: $!";
my $header = <FILE>;
while(<FILE>){
	chomp $_;
	my @cols = split(/\t/,$_);
	if ((!(defined $best_ids{$cols[0]}))or($best_ids{$cols[0]} < $cols[2])){
		$best_ids{$cols[0]} = $cols[2];
	}
}
close FILE or die "Error closing $file: $!";

# compile an array of the data
my @data;
foreach my $id (keys %best_ids){
	push(@data,$best_ids{$id});
}

my $graph = new GD::Graph::histogram(400,600);
$graph->set( 
	x_label         => 'Percent Identity of best match',
	y_label         => 'Frequency',
	title           => 'Histogram of wateredBLAST results for: '.$file,
	x_labels_vertical => 1,
	bar_spacing     => 0,
	shadow_depth    => 1,
	shadowclr       => 'dred',
	transparent     => 0,
	histogram_bins	=> 25,
	histogram_type	=> 'percentage'
	
) or warn $graph->error;

# plot
my $gd = $graph->plot(\@data) or die $graph->error;

# create the image
open(IMG, ">$output") or die $!;
binmode IMG;
print IMG $gd->png;

