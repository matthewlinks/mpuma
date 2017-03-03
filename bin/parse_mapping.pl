#!/usr/bin/perl -w
use strict;

use Getopt::Long;

my ($mappingDir,$mappingFile,$assemblyDir,$VERBOSE);

my @options = (
	'd=s',	\$mappingDir,
	'f=s',	\$mappingFile,
	'a=s',	\$assemblyDir,
	'v=s',	\$VERBOSE

);
&GetOptions(@options);
die "Usage: $0 -d mappingDir sff1 sff2 ... sffN" unless(
	defined $mappingDir
);
my @sffFiles = @ARGV; 
die "Error need at least 1 sff file" unless (@sffFiles > 0);
$mappingFile = join('/',$mappingDir,'454ReadStatus.txt') unless defined $mappingFile;
$assemblyDir = join('/',$mappingDir,'Sub-assemblies') unless defined $assemblyDir;
mkdir($assemblyDir) or die "Error creating $assemblyDir" unless (-e $assemblyDir);


my $mapping = parse_mapping($mappingFile);

foreach my $OTUbin (sort {return @{$mapping->{$b}} <=> @{$mapping->{$a}}} keys %{$mapping}){
	my $count = @{$mapping->{$OTUbin}};
	warn join("\t",$OTUbin,$count) if $VERBOSE;
	my $aDir = join('/',$assemblyDir,$OTUbin);
	if (!(-e $aDir)){
		mkdir($aDir) or die "Error creating directory $aDir: $!";
	}

	# make sure that there is a read file in place
	my $readFile = join('/',$aDir,'reads.ids');
	open(FILE,">$readFile") or die "Error opening $readFile: $!";
	foreach my $id (@{$mapping->{$OTUbin}}){
		print FILE $id,"\n";
	}
	close FILE or die "Error closing $readFile: $!";

	# create an SFF from those reads
	my $outSFF = join('/',$aDir,'reads.sff');
	system('sfffile',
		'-o',	$outSFF,
		'-i',	$readFile,
		@sffFiles) == 0 or die "Error creating $outSFF: $!";
}

exit 0;

sub parse_mapping {
	my $file = shift;
	die "Error lost file" unless defined $file;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my %map;
	my %seenReads;
	while(<FILE>){
		chomp $_;
		my @parts = split(/\t/,$_);
		if (($parts[1] eq 'Partial')or($parts[1] eq 'Full')){
			next if defined $seenReads{$parts[0]}->{$parts[4]};
			push(@{$map{$parts[4]}},$parts[0]);
			$seenReads{$parts[0]}->{$parts[4]}++;
		}elsif(($parts[1] eq 'Repeat')or($parts[1] eq 'Unmapped')){
			next if defined $seenReads{$parts[0]}->{'Unknown'};
			push(@{$map{'Unknown'}},$parts[0]);
			$seenReads{$parts[0]}->{'Unknown'}++;
		} # else Chimeric, or TooShort
	}
	close FILE or die "Error closing $file: $!";
	return \%map;
}
