#!/usr/bin/perl -w
use strict;
use lib "$ENV{mPUMA}/lib";
use mPUMA;
use Getopt::Long;

my ($workingDir,$libDir);

my @options = (
	'w=s',	\$workingDir,
	'l=s',	\$libDir,
);
&GetOptions(@options);

die "Usage: $0 -l libDir -w workingDir"
unless (defined $libDir and defined $workingDir);

my $libs = get_subDirs($libDir);
my @divFiles;
foreach my $lib (@{$libs}){
	my $diversityFile = join('/',$libDir,$lib,'diversity');
	push(@divFiles,$diversityFile);
}

    
# Create a list of isotigs that exist in the diversity file
my $OTU_nuc_list_file = mPUMA::generateIsotigList( $workingDir,\@divFiles);

# Creates template technology files for OTU extraction used by the convertTech sub routine.
my $technology_nuc_file = mPUMA::generateTechFile( $OTU_nuc_list_file, 10000);

# Converts the template files into technology (*.tsv) files usable by gene_spring.
# Checks the convertedTech output formated file for errors.
my ($genespring_nuc_dir,$status) = mPUMA::convertTech( $workingDir, $technology_nuc_file,$libDir,$libs);




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
