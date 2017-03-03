#! /usr/bin/perl -w
use strict;
use Getopt::Long;
my ($divDir, $output);
my @options = (
	'd=s',	\$divDir,
	'o=s',	\$output,
);
&GetOptions(@options);
usage() unless (
	defined $divDir
	and defined $output
);

sub usage {

die <<"USAGE";

	Usage: $0 -d genespring_dir -o outputUnifrac

	   *   genespring_dir the genespring diversity file directory. 
	   *   outputUnifrac name of the unifrac input file.

USAGE
}
my ($DIR, @divListing, $unifracFile, @isotigList, %seen, @numIsotigs);
@divListing = ();

opendir($DIR, $divDir);
while(my $divfiles = readdir($DIR)){
	if ($divfiles =~ m/.+\.txt/ and $divfiles !~ m/.+\.txt~/ ){
		push(@divListing, $divfiles);
	}
}
closedir($DIR);

@isotigList = ();
@numIsotigs = ();

my %unifrac;
foreach my $divs (@divListing){
	my $libName;
	if($divs =~ m/(.+)\.txt/){
		$libName = $1;
	}else{
		die "Error with filename";
	}

	my $genespring_file = join('/', $divDir, $divs);
	open(INFILE, "<", $genespring_file) or die $!; 
	while(<INFILE>){
# 		isotig02561	0	0	0	something
		if($_ =~ m/^(isotig\d+)\t(\d+|\d+.\d+)\t(\d+|\d+.\d+)\t(\d+|\d+.\d+)\t.+\n/){
			my $isotig = $1;
			push(@isotigList, $isotig);
			my $scaledNum = $3;
			push(@{$unifrac{$libName}}, $scaledNum);
		}
	}
	close INFILE or die $!;
}
@isotigList = grep { ! $seen{$_} ++ } @isotigList;
@isotigList = sort(@isotigList);

$unifracFile = join('/', $output, 'unifrac_freq_table.txt');

unless(-s $unifracFile){
	open OUTFILE, ">", $unifracFile or die $!;
	my $i = 0;
	foreach my $lib (keys %unifrac){
		$i = 0;
		foreach my $isotig (@isotigList){

			print OUTFILE join("\t", $isotig, $lib, ${$unifrac{$lib}}[$i]),"\n";
			$i++;
		}
	}
	close OUTFILE or die $!;
}

exit 0;
