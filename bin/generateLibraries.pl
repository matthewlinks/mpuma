#!/usr/bin/perl -w
use strict;

use Getopt::Long;
my $Dir;
my @options = (
	'd=s',	\$Dir,	
);
&GetOptions(@options);

my @SFF_List = generateSFFLibs($Dir, $Dir . "sff");

exit 0;

sub generateSFFLibs{
	my ($storageDir, $libraryList, $DIR, $SFFsDir, @listingSFFs, @SFF_List);
	$storageDir = shift;
        die "Error lost the SFF storage directory" unless defined $storageDir;
	$libraryList = $storageDir . "library-sff-list";
	open OUTFILE, ">", $libraryList or die $!;
	foreach my $libraryDir (@_){
		@listingSFFs = ();
		opendir($DIR, $libraryDir);
		while(my $SFFfiles = readdir($DIR)){
			if ($SFFfiles =~ m/.+MID\d+.sff/){
				push(@listingSFFs, $SFFfiles);
			}
		}
		closedir($DIR);
	
		foreach my $SFFs (@listingSFFs){
			my $SFFsPath = $libraryDir . "/" . $SFFs;
			$SFFs =~ s/.sff//g;
			my $line = "$SFFs\t$SFFsPath\n";
			print OUTFILE $line;
			print "$SFFsPath\n";
			push(@SFF_List, $SFFsPath);
		}
	}
	close OUTFILE or die $!;
	return @SFF_List;
}
