#!/usr/bin/perl -w
$| = 1;
use strict;
use JSON;
use Getopt::Long;
my ($inputFile,$outputFile);
# what are the command line options
my @options = (
	'i=s',	\$inputFile,
	'o=s',	\$outputFile,
);
&GetOptions(@options);

# we want to ensure that there is a file as input
# for our output we will print to a file when its there
# but its not required that we have an output file
# if the output file isn't there then we'll just print
# to STDOUT
die "Usage: $0 -i input -o output" 
unless (defined $inputFile);

# read the file into a variable
# N.B. this variable contains JSON and not a PERL object yet
my $json = read_biom($inputFile);

# JSON => PERL hash reference
my $object = decode_json($json);
# now $object is just a PERL hash that we can play around with 

# before we play around with the object let's setup our output stream so we can just print and not care anymore about where its going
# this line is going to redefine STDOUT if we have a filename to try and write to
open(STDOUT,">$outputFile") or die "Error opening $outputFile for writting: $!" if (defined $outputFile);

# In this section do something with $object


# Example - check for an OTUtree (should be newick) and if there is one output it
# you may or maynot want to add a newline after this...
#if (defined $object->{OTUtree}){
#	print $object->{OTUtree};
#}

# Example - produce a FASTA file of the OTU sequences
## print the sequence pretty - $WIDTH chars wide
#my $WIDTH = 60;
#foreach my $otu (@{$object->{rows}}){
#	print ">",$otu->{id},"\n";
## 	If you want the sequence to all be on one line then just go 
##	print $otu->{metadata}->{sequence},"\n";
##	otherwise you can print it with a max width (e.g. 60)
#	my $seq = "\n";
#	if (defined $otu->{metadata}->{sequence}){
#		$seq = $otu->{metadata}->{sequence};
#		if (length($seq) > 0){
#			$seq =~ s/(.{1,$WIDTH})/$1\n/g;	
#		}
#	}
#	print $seq;
#}

# Example - produce a tab-delimited OTU table
# we need to build up the order of the OTUs because the data table doesn't record the OTU IDs in situ
my @OTUorder; 
foreach my $otu (@{$object->{rows}}){
	push(@OTUorder,$otu->{id});
}

# we need to build up the order of the Libraries (samples) because the data table doesn't record the OTU IDs in situ
my @LibOrder;
foreach my $lib (@{$object->{columns}}){
	push(@LibOrder,$lib->{id});
}

my $DELIM = "\t";
# write a header line 
print join($DELIM,'',@LibOrder),"\n";

# print out a line for the abundance of each OTU (row) by library (cols)
for(my $otuIndex = 0; $otuIndex < @OTUorder; $otuIndex++){
	my $otu = $OTUorder[$otuIndex]; # the name of the OTU
	die "Error with the abundance for $otu" unless defined $object->{data}->[$otuIndex];
	my $abundances = $object->{data}->[$otuIndex]; 
	print join($DELIM,$otu,@{$abundances}),"\n";
}

# we are done so we can just exit
exit 0;


sub read_biom {
	my $file = shift;
	die "Error lost file" unless defined $file;
	open(my $fh, "<", $file) or die "Error opening $file for reading: $!";
	my $j = <$fh>;
	return $j;
}



