#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($read_file,$num_reads,$output);
my @options = (
	'r=s',	\$read_file,
	'n=s',	\$num_reads,
	'o=s',	\$output
);
&GetOptions(@options);
usage() unless(defined $read_file and defined $num_reads);

# redirect to a file if requested
if (defined $output){
	open(STDOUT,">$output") or die "Error opening $output for writting: $!";
}

my $accs = read_accs($read_file);
while($num_reads > 0){
	my $acc = pick_at_random($accs);
	print $acc,"\n";
	$num_reads--;
	my $numRemaining = @{$accs};
	last if ($numRemaining < 1);
}
exit;

sub usage{
	die "Usage: $0 -r read_file -n number";
}

sub pick_at_random {
	my $hat = shift;
	die "Error lost names" unless defined $hat;

	# make sure there is a name in the hat otherwise return undef
	my $numInHat = @{$hat};
	return undef if ($numInHat < 1);

	# pick a name from the hat
	my $name;
	my $random_int;
	while(!defined $name){
		$random_int = int(rand(@{$hat}));	
		$name = $hat->[$random_int];
	}

	# remove that name from the hat
	delete $hat->[$random_int];
	
	# return the name from the hat
	return $name;
}

sub read_accs {
	my $file = shift;
	die "Error lost filename" unless defined $file;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my %reads;
	while(<FILE>){
		chomp $_;
		my ($acc) = split(/\s+/,$_);
		$reads{$acc}++;
	}
	close FILE or die "Error closing $file: $!";
	
	my @total = keys %reads;
	return \@total;
}
