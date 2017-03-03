#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
my ($dir,$extension);

my @options = (
	'd=s',	\$dir,
	'e=s',	\$extension
);
&GetOptions(@options);
$extension = 'phylum' unless defined $extension;
die "Usage: $0 -d dir -e ext" unless (
	defined $dir
	and defined $extension
);


my %results;
my %taxa;
my @libs;

opendir(DIR,$dir) or die "Error opening $dir: $!";
while(my $file = readdir(DIR)){
	next if ($file =~ /^\./);
	next unless ($file =~ /\.$extension$/);
	my $lib = $file;
	if ($file =~ /^(.*)\.$extension$/){
		$lib = $1;
	}
	push(@libs,$lib);
	open(FILE,"<$dir/$file") or die "Error opening $dir/$file: $!";
	my $header = <FILE>;
	while(<FILE>){
		chomp $_;
		my ($t,$percent) = split(/\t/,$_);
		$taxa{$t}++;
		$results{$lib}->{$t} = $percent;
	}
	close FILE or die "Error closing $file: $!";
}
closedir DIR or die "Error closing $dir: $!";

my @orderedTaxa = sort {return $a cmp $b} keys %taxa;
print join("\t",'Library',@orderedTaxa),"\n";
foreach my $lib (@libs){
	my @values;
	foreach my $t (@orderedTaxa){
		my $val = 0;
		$val = $results{$lib}->{$t} if defined $results{$lib}->{$t};
		push(@values,$val);
	}
	print join("\t",$lib,@values),"\n";
}

exit 0;
















