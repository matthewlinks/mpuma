#!/usr/bin/perl -w
use strict;
use lib "$ENV{mPUMA}/lib";
use mPUMA;
use Getopt::Long;

my ($workingDir,$libDir);

my @options = (
	'w=s',	\$workingDir,
	'l=s',	\$libDir
);
&GetOptions(@options);

die "Usage: $0 -l libDir -w workingDir"
unless (defined $libDir and defined $workingDir);

my $libs = get_subDirs($libDir);

mPUMA::generateInputMOTHUR($workingDir,$libDir,$libs);
my ($mothur_nuc_file, $mothur_nuc_files_dir) = mPUMA::generateInputMOTHUR($workingDir,$libDir,$libs);
mPUMA::processInputMOTHUR( $mothur_nuc_file, $mothur_nuc_files_dir);

sub get_subDirs{
	my $dir = shift;
	opendir(DIR,$dir) or die "Error opening $dir: $!";
	my @dirs;
	while(my $file = readdir(DIR)){
		next if ($file =~ /^\./);
		next unless (-d "$dir/$file");
		push(@dirs,$file);
	}
	closedir(DIR);
	return \@dirs;
}
