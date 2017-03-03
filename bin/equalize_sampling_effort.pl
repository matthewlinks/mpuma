#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#use File::Basname;

my ($libraryDir,$outputDir);
my @options = (
	'l=s',	\$libraryDir,
	'o=s',	\$outputDir
);
&GetOptions(@options);

die "Usage: $0 -l libraryDir -o outputDir"
unless (
	defined $libraryDir
	and defined $outputDir
);

# figure out what libraries exist based on reads.ids files in those sub directories
my $libs = get_libraries($libraryDir);

# figure out what the minimum number of reads is
my $minNumberOfReads = find_minimum_number_of_reads($libs);

# make sure the output directory exists
if (!(-e $outputDir)){
	mkdir ($outputDir) or die "Error creating directory $outputDir: $!";
}

# downsample
my @downSampledFiles;
foreach my $lib (keys %{$libs}){
	my $readFile = $libs->{$lib};

	# down Sample the library
	my $downSampledFile = join('/',$outputDir,$lib . '.downSampledIDs');
	die "Error $downSampledFile already exists so aborting" if (-e $downSampledFile);
	system("$ENV{mPUMA}/bin/pick_random_reads.pl",
		'-n',	$minNumberOfReads,
		'-r',	$readFile,
		'-o',	$downSampledFile) == 0 or die "Error picking random reads from $readFile: $!";

	push(@downSampledFiles,$downSampledFile);
}

print join(' ',@downSampledFiles),"\n";

exit 0;


sub get_libraries {
	my $dir = shift;
	die "Error lost dir" unless defined $dir;
	opendir(DIR,$dir) or die "Error opening $dir for reading: $!";
	my %libraries;
	while(my $lib = readdir(DIR)){
		my $readFile = join('/',$dir,$lib,'reads.ids');
		$libraries{$lib} = $readFile if (-e $readFile);
	}
	closedir(DIR) or die "Error closing $dir: $!";
	return \%libraries;
}

sub find_minimum_number_of_reads {
	my $libs = shift;
	die "Error lost library hash" unless defined $libs;
	my $min;
	foreach my $lib (keys %{$libs}){
		my $readFile = $libs->{$lib};
		my $num = 0;
		open(FILE,"<$readFile") or die "Error opening $readFile for reading: $!";
		$num++ while(<FILE>);
		close FILE or die "Error closing $readFile: $!";
		if (! defined $min){
			$min = $num;
		}elsif($num < $min){
			$min = $num;
		}
	}
	return $min;
}
