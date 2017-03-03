#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $working_dir;
my @options = (
	'd=s',	\$working_dir,
);
&GetOptions(@options);

usage() unless (
	defined $working_dir
);

sub usage {
	die "Usage: $0 -d $working_dir";
}

my ($DIR, @divListing, $megan_dir);

@divListing = ();

$megan_dir = join('/', $working_dir);
unless(-d $megan_dir){
		mkdir($megan_dir, 0777) or die "Can't make directory: $!";
}

opendir($DIR, $megan_dir);
while(my $megan_files = readdir($DIR)){
	if ($megan_files =~ m/(.+)\.megan.blastx$/){
		my $lib = $1;
		my $rmalib = $lib . '.rma'; 
		my $rmafile = join('/', $working_dir, $rmalib);
		unless(-e $rmafile){
			my $meganCmd = "MEGAN +g -x \"set dir=$megan_dir; import blastfile=$lib.megan.blastx fastafile=$lib.megan.fasta meganfile=$lib.rma minscore=35.0; quit;\"";
			system(
			$meganCmd,
			) == 0 or die "Error calling $meganCmd: $!";
		}
		else{
			print "$rmafile has been previously created....\n";
		}
	}
}
closedir($DIR);

exit 0;
