#!/usr/bin/perl -w
use strict;
use File::Basename;
while(<>){
	chomp $_;
	my @parts = split(/\//,$_);
	my ($pool,$mid);
	for(my $i = 0; $i < @parts; $i++){
		if ($parts[$i] =~ /^([BA][CD]\-\d+)\-sffs$/){
			$pool = $1;
			if ($parts[$i+1] =~ /^\w+\-(MID\d+)[\-\.]/){
				$mid = $1;
				last;
			}else{
				die "Error can parse pool but not mid";
			}
		}
	}
	die "Error parsing pool and mid from $_" unless (defined $pool and defined $mid);
	my $sff = basename($_);
	print join("\t",join('-',$pool,$mid),$sff),"\n";

}
exit 0;
