#!/usr/bin/perl -w
use Statistics::Descriptive;
use strict;
use Getopt::Long;
my ($dir,$output,$col,$selectOTUfile,$header);
my @options = (
	'd=s',	\$dir,
	'o=s',	\$output, 
	'c=s',	\$col,
	'i=s',	\$selectOTUfile,
	'h=s',	\$header
);
&GetOptions(@options);
$col = 1 unless defined $col;
$header = 0 unless defined $header;
die "Usage: $0 -d dir -o output -c 1" unless defined $dir;

my %resultsByOTU;
my %libs;

my $selectResults = read_ids($selectOTUfile);

opendir(DIR,$dir) or die "Error opening $dir: $!";
while(my $file = readdir(DIR)){
	next if ($file =~ /^\./);
	open(FILE,"<$dir/$file") or die "Error opening $file: $!";
	my ($lib,$rep);
	if ($file =~ /^(.*)\_(\d+)/){
		$lib = $1;
		$rep = $2;
	}else{
		die "Error cannot figure out the replicate number from $file";
	}
	$libs{$lib}++;
	my $header = <FILE> if ($header);
	while(<FILE>){
		chomp $_;
		my ($otu,@cols) = split(/\t/,$_);
		if(defined $selectOTUfile){
			if (defined $selectResults->{$otu}){
				die "Error $otu is already defined for $lib $rep" if defined $resultsByOTU{$otu}->{$lib}->{$rep};
				$resultsByOTU{$otu}->{$lib}->{$rep} = $cols[$col - 1];
			} # else ignore it
		}else{
			# keep everything
				die "Error $otu is already defined for $lib $rep" if defined $resultsByOTU{$otu}->{$lib}->{$rep};
				$resultsByOTU{$otu}->{$lib}->{$rep} = $cols[$col-1];
		}
	}
	close FILE or die "Error closing $file";
}

my @orderedLibs = sort {return $a cmp $b} keys %libs;

# write to a file if specified
if (defined $output){
	open(STDOUT,">$output") or die "Error opening $output: $!";
}

print join("\t",'OTU',@orderedLibs),"\n";

foreach my $otu (sort {return $a cmp $b} keys %resultsByOTU){
	my @medians;
	foreach my $lib (@orderedLibs){
		my @values;
		foreach my $rep (keys %{$resultsByOTU{$otu}->{$lib}}){
			push(@values,$resultsByOTU{$otu}->{$lib}->{$rep});
		}
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@values);
		my $median = $stat->median();
		push(@medians,$median);
	}
	print join("\t",$otu,@medians),"\n";
}

exit 0;


sub read_ids {
	my $file = shift;
	return undef unless defined $file;
	my %ids;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	while(<FILE>){
		chomp $_;
		my ($id) = split(/\s+/,$_);
		$ids{$id}++;
	}
	close FILE or die "Error closing $file: $!";
	return \%ids;
}

