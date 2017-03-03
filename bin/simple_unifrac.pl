#!/usr/bin/perl -w
use strict;
use lib "$ENV{mPUMA}/lib";
use mPUMA;
use Getopt::Long;

my ($workingDir,$libDir,$orientedFile);

my @options = (
	'w=s',	\$workingDir,
	'l=s',	\$libDir,
	'o=s',	\$orientedFile,
);
&GetOptions(@options);

die "Usage: $0 -l libDir -w workingDir -o orientedFile"
unless (defined $libDir and defined $workingDir and defined $orientedFile);

my $libs = get_subDirs($libDir);
my @divFiles;

foreach my $lib (@{$libs}){
	my $diversityFile = join('/',$libDir,$lib,'diversity');
	push(@divFiles,$diversityFile);
}

my $nucUnifracDir = mPUMA::generateInputUNIFRAC( $workingDir, \@divFiles, $orientedFile,'Nucleotide');

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
