#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($dir,$output,$suffix,$min,$exclude);
my @options = (
	'd=s',	\$dir,
	'o=s',	\$output,
	's=s',	\$suffix,
	'm=s',	\$min,
	'e=s',	\$exclude,
);
&GetOptions(@options);

usage() unless(
	defined $dir
	and defined $output
);
$suffix = 'txt' unless defined $suffix;
$min = 0 unless defined $min;
$exclude = 0 unless defined $exclude;

sub usage {
	die "Usage: $0 -d dir -o output -s 'txt'";
}

my $files = find_files($dir,$suffix);

# two hashes to store the abundance and the OTUnames we find in the files
my %OTUabundance;
my %OTUnames;

# read in the data from $file
foreach my $file (@{$files}){
	my $path = join('/',$dir,$file);
	open(FILE,"<$path") or die "Error opening $path for reading: $!";
	while(<FILE>){
		chomp $_;
		my @parts = split(/\s+/,$_);
		my $id = $parts[0];
		next if (($exclude)and($id !~ /^isotig/));

		my $actualCount = $parts[3];
		next unless $actualCount >= $min; # allow the user to specify a minimum abundance level to exclude things like singletons as per their preference
		$OTUabundance{$file}->{$id} = $actualCount;
		$OTUnames{$id}++;
	}
	close FILE or die "Error closing $file: $!";
}

# write out the output
open(OUTPUT,">$output") or die "Error opening $output for writting: $!";
my @sortedOTUnames = sort {return $a cmp $b} keys %OTUnames; # sort the OTUnames so that we can ensure that the columns in the MOTHUR input file will line up accordingly
my $numberOfOTU = @sortedOTUnames;

# I wish mothur would allow some type of commenting in the file format so there there could be some reuse of this datafile...
#print OUTPUT "#",join("\t",'NULL','LIBRARY NAME','#OTUs',@sortedOTUnames),"\n";
foreach my $file (keys %OTUabundance){
	if ($file =~ /^(.*)\.$suffix$/){
		my $lib = $1;

		# create the ordered set of abundance values
		my @orderedAbundances;
		foreach my $id (@sortedOTUnames){
			my $abundance = 0;	# if there was no abundance (or it was below $min) then set this to 0
			$abundance = $OTUabundance{$file}->{$id} if defined $OTUabundance{$file}->{$id};
			push(@orderedAbundances,$abundance);
		}

		# ensure that there are the same number of abundances as we expect!!!
		my $numAbundances = @orderedAbundances;
		die "Error $numAbundances != $numberOfOTU and so this would corrupt MOTHUR input file"
			unless ($numAbundances == $numberOfOTU);

		print OUTPUT join("\t",'NA',$lib,$numberOfOTU,@orderedAbundances),"\n";
	}else{
		die "Error $file does not seem to have the suffix $suffix";
	}
}

close OUTPUT or die "Error closing $output: $!";

exit 0;


sub find_files {
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	my $suffix = shift;
	die "Error lost pattern" unless defined $suffix;

	opendir(DIR,$dir) or die "Error opening $dir: $!";

	my @files;
	while(my $file = readdir(DIR)){
		push(@files,$file) if ($file =~ /\.$suffix$/);
	}
	closedir(DIR) or die "Error closing $dir: $!";
	return \@files;
}

