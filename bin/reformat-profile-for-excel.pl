#!/usr/bin/perl -w
use strict;
use File::Basename;

my %terms;
my %results;
foreach my $file (@ARGV){
	my $basename = basename($file);
	if ($basename =~ /^(.*)\..*$/){
		$basename = $1;
	}
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $header = <FILE>;
	while(<FILE>){
		chomp $_;
		my ($term,$val) = split(/\t/,$_);
		$results{$basename}->{$term} = $val;
		$terms{$term}++;
	}
	close FILE or die "Error close $file: $!";
}
my @termOrder = sort { return $a cmp $b } keys %terms;

print join("\t",'File',@termOrder),"\n";

foreach my $file (sort {return $a cmp $b} keys %results){
	my @vals;
	foreach my $term (@termOrder){
		my $v = 0;
		$v = $results{$file}->{$term} if (defined $results{$file}->{$term});
		push(@vals,$v);
	}
	print join("\t",$file,@vals),"\n";
}
exit 0;
