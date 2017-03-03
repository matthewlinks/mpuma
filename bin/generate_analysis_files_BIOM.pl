#!/usr/bin/perl -w
use strict;
use POSIX qw(strftime);

use File::Basename;
use Getopt::Long;

my ($otuFile,$diversityDir,$suffix,$col,$header,$outFile,$wateredBLASTfile,$wateredBLASTheader);
my @options = (
	'i=s',	\$otuFile,
	'd=s',	\$diversityDir,
	's=s',	\$suffix,
	'h=s',	\$header,
	'c=s',	\$col,
	'o=s',	\$outFile,
	'w=s',	\$wateredBLASTfile,
	'wh=s',	\$wateredBLASTheader
);
&GetOptions(@options);
$suffix = '.txt' unless defined $suffix;
$header = 1 unless defined $header;
$wateredBLASTheader = 0 unless defined $wateredBLASTheader;
$col = 3 unless defined $col;

my $files = find_files($diversityDir,$suffix);
my %data;
my $otus = read_otu_file($otuFile);
my $matches = read_wateredBLAST($wateredBLASTfile,$wateredBLASTheader) if defined $wateredBLASTfile;

foreach my $lib (keys %{$files}){
	my $file = $files->{$lib};
	open(FILE,"<$file")or die "Error opening $file for reading: $!";
	my $header = <FILE> if ($header);
	while(<FILE>){
		my ($OTU,@cols) = split(/\t/,$_);
		if (defined $otuFile){
			next unless defined $otus->{$OTU};
		}
		my $val = $cols[$col - 1]; # columns are 0-indexed and we loose 1 to the otu label
		die "Error multiple values for $OTU in $lib" if defined ($data{$lib}->{$OTU});
		$data{$lib}->{$OTU} = $val;
#		$otus->{$OTU}++;
	}
	close FILE or die "Error closing $file";
}

my @orderedLibs = sort{return $a cmp $b} keys %{$files};
my @orderedOTUs = sort{return $a cmp $b} keys %{$otus};
my $num = @orderedOTUs;

sub read_wateredBLAST{
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $header = shift; 
	$header = 0 unless defined $header;
	
	my %matches; 
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $h = <FILE> if ($header);

	while(<FILE>){
		chomp $_;
		my ($query,$qDesc,$percentMatch,$lengthMatch,$subject,$strand) = split(/\t/,$_);
		next if defined $matches{$query}; # only report 1st match
		$matches{$query} = join(',',$subject,$percentMatch,$lengthMatch,$strand);
	}
	close FILE or die "Error closing $file: $!";
	return \%matches;
}
sub get_date{
	my $now = time();
	return strftime("%Y-%m-%dT%H:%M:%S", localtime($now));
}

open(STDOUT,">$outFile") if defined $outFile;

my $numOTUs = @orderedOTUs;
my $numLibs = @orderedLibs;
my $date = get_date();

# initial block
print <<EOF;
{
	"id":null,
	"format": "Biological Observation Matrix 0.9.1-dev",
	"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
	"type": "OTU table",
	"generated_by": "mPUMA revision 1",
	"date": "$date",
	"rows":[
EOF

# print the row information for each OTU
my @rowElements;
for(my $i =  0; $i < $numOTUs;$i++){
	my $otu = $orderedOTUs[$i];
	if (defined $wateredBLASTfile){
		my $match = 'null';
		$match = $matches->{$otu} if (defined $matches->{$otu});
		print "\t\t{\"id\":\"$otu\", \"metadata\":\"$match\"}";
	}else{
		print "\t\t{\"id\":\"$otu\", \"metadata\":\"null\"}";
	}
	if ($i < ($numOTUs -1)){
		print ",";
	}
	print "\n";
}
print join(",\n",@rowElements),"\n";

print <<EOF;
	],
	"columns": [
EOF
# print the library information 
my @libElements;   
for(my $i =  0; $i < $numLibs;$i++){
	my $lib = $orderedLibs[$i];
	print "\t\t{\"id\":\"$lib\", \"metadata\":null}";
	if($i < ($numLibs - 1)){
		print ",";
	}
	print "\n";
}

print <<EOF;
	],
	"matrix_type": "dense",
	"matrix_element_type": "float",
	"shape": [$numOTUs,$numLibs],
	"data":	[
EOF

# print out the data 
for(my $i = 0; $i < $numOTUs;$i++){
	my $otu = $orderedOTUs[$i];
	my @vals;
	foreach my $lib (@orderedLibs){ 
		my $v = 0;
		$v = $data{$lib}->{$otu} if defined $data{$lib}->{$otu};
		push(@vals,$v);
	}
	print "\t\t[",join(',',@vals),"]";
	if ($i < ($numOTUs -1)){
		print ",";
	}
	print "\n";
}
print <<EOF;
		]
}
EOF

exit 0;

sub read_otu_file {
	my $file = shift;
	die "Error lost file" unless defined $file;
	open(FILE,"<$file");
	my %OTUs;
	while(<FILE>){
		chomp $_;
		my ($id) = split(/\s+/,$_);
		my $alias = $id;
		if ($id =~ /(.*)\(.*\)$/){
			$id = $1;
			# there is something in brackets which we would like to see in the resulting images but it will mess up the accounting if its used to match OTU labels
		}	
		$OTUs{$id} = $alias; 
#		$OTUs{$id} = 1;
	}
	close FILE or die "Error closing $file: $!";
	return \%OTUs;
}


sub find_files{
	my $dir = shift;
	die "Error lost dir" unless defined $dir;
	my $suffix = shift;
	die "Error lost suffix" unless defined $suffix;

	my %files;
	opendir(DIR,$dir) or die "Error opening $dir: $!";
	while(my $file = readdir(DIR)){
		next if ($file =~ /^\./);
		warn "Checkin $file vs $suffix";
		if ($file =~ /(.*)$suffix$/){
			my $lib = $1;
			warn "$lib => $dir/$file";
			$files{$lib} = join('/',$dir,$file);
		}		
	}
	close DIR;
	return \%files;
}
