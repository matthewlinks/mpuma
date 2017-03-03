#!/usr/bin/perl -w
use strict;
use lib "$ENV{mPUMA}/lib";
use mPUMA;
use Getopt::Long;

my ($workingDir,$libDir,$TAXONOMICTERM,$classifierFile);

my @options = (
	'w=s',	\$workingDir,
	'l=s',	\$libDir,
	't=s',	\$TAXONOMICTERM,
	'c=s',	\$classifierFile
);
&GetOptions(@options);

die "Usage: $0 -l libDir -w workingDir"
unless (defined $libDir and defined $workingDir and defined $classifierFile);
$TAXONOMICTERM = 'phylum' unless defined $TAXONOMICTERM;

my $libs = get_subDirs($libDir);

foreach my $lib (@{$libs}){
	my $diversityFile = join('/',$libDir,$lib,'diversity');
	my $meganFile = mPUMA::generateClassifierInputMEGAN( $workingDir, $lib, $diversityFile, $classifierFile);
	mPUMA::generateClassifierProfilesFromMegan($workingDir, $lib,$meganFile,$TAXONOMICTERM) if (defined $TAXONOMICTERM);
}

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
