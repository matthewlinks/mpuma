#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($libDir,$readFile,$subSampleDir,$min,$VERBOSE);

my @options = (
	'l=s',	\$libDir,
	'r=s',	\$readFile,
	's=s',	\$subSampleDir,
	'm=s',	\$min,
	'v=s',	\$VERBOSE
);
&GetOptions(@options);

die "Usage: $0 -l libDir -r reads.ids -s outputDir -m minimum"
unless (
	defined $libDir
	and defined $subSampleDir
);
$readFile = 'reads.ids' unless defined $readFile;
$VERBOSE = 0 unless defined $VERBOSE;

my $lib2reads = find_libraries($libDir);

# find the minimum value
if (!defined $min){
	warn "Calculating minimum # of reads in the current ID lists" if $VERBOSE;
	foreach my $lib (keys %{$lib2reads}){
		if (!defined $min){
			$min = $lib2reads->{$lib};
		}else{
			$min = $lib2reads->{$lib} if ($min > $lib2reads->{$lib});
		}
	}
	warn "The minimum to sample down to was $min" if $VERBOSE;
} # else use whatever we were told ===> N.B. this may explode or have unexpected consequences if you ask for more reads than are in a library to begin with (because random sampling is done WITHOUT replacement)

# make sure there is a directory for the subsampling exercise
ensure_directory($subSampleDir);

# actually sub sample by read ID
foreach my $lib (keys %{$lib2reads}){
	warn "Generating subsampling data for $lib" if $VERBOSE;
	my $src = join('/',$libDir,$lib,$readFile);
	my $dest = join('/',$subSampleDir,join('.',$lib,'subsampled-to',$min,'txt'));
	system("$ENV{mPUMA}/bin/pick_random_reads.pl",
		'-n',	$min,
		'-r',	$src,
		'-o',	$dest) == 0 or die "Error sub-sampling $src to $min reads and writting to $dest: $!";
}

exit 0;

sub ensure_directory {
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	if (-e $dir){
		die "Error $dir exists but is NOT a directory" unless (-d $dir);
	}else{
		mkdir($dir) or die "Error creating directory $dir";
	}
}

sub find_libraries {
	my $dir = shift;
	die "Error lost dir" unless defined $dir;
	warn "Looking for libraries in $dir" if $VERBOSE;
	opendir(DIR,$dir) or die "Error opening $dir: $!";
	my %libCounts;
	while(my $file = readdir(DIR)){
		warn "Checking ",join('/',$dir,$file,$readFile) if $VERBOSE;
		if (-s join('/',$dir,$file,$readFile)){
			$libCounts{$file} = count_lines(join('/',$dir,$file,$readFile));
		}
	}
	closedir DIR or die "Error closing $dir: $!";
	return \%libCounts;
}

sub count_lines {
	my $file = shift;
	die "Error lost file" unless defined $file;
	open(FILE,"<$file") or die "Error opening $file: $!";
	my $n = 0;
	$n++ while(<FILE>);
	close FILE or die "Error closing $file: $!";
	return $n;
}

