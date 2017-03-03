#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

my ($otuFile,$diversityDir,$suffix,$col,$header,$outFile,$transpose);
my @options = (
	'i=s',	\$otuFile,
	'd=s',	\$diversityDir,
	's=s',	\$suffix,
	'h=s',	\$header,
	'c=s',	\$col,
	'o=s',	\$outFile,
	't=s',	\$transpose
);
&GetOptions(@options);
$suffix = '.txt' unless defined $suffix;
$header = 1 unless defined $header;
$col = 3 unless defined $col;
$transpose = 0 unless defined $transpose;

my $files = find_files($diversityDir,$suffix);
my %data;
my $otus = read_otu_file($otuFile);

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
#		warn "$lib => $OTU = $val";
#		$otus->{$OTU}++;
	}
	close FILE or die "Error closing $file";
}

my @orderedLibs = sort{return $a cmp $b} keys %{$files};
my @orderedOTUs = sort{return $a cmp $b} keys %{$otus};
my $num = @orderedOTUs;

open(STDOUT,">$outFile") if defined $outFile;
if ($transpose){
	print join("\t",@orderedLibs),"\n";
	foreach my $otu (@orderedOTUs){
		my @elements;				
		push(@elements,$otus->{$otu});
			
		foreach my $lib (@orderedLibs){
			my $val = 0;
			$val = $data{$lib}->{$otu} if defined $data{$lib}->{$otu};
			push(@elements,$val);
		}
		print join("\t",@elements),"\n";
	}
}else{
	my @aliasedOTUs;
	foreach my $otu (@orderedOTUs){
		push(@aliasedOTUs,$otus->{$otu});
	}
	print join("\t",@aliasedOTUs),"\n";
	foreach my $lib (@orderedLibs){
		my @elements;
		push(@elements,$lib);
		
		foreach my $otu (@orderedOTUs){
			my $val = 0;
			$val = $data{$lib}->{$otu} if defined $data{$lib}->{$otu};
			push(@elements,$val);
		}		
		print join("\t",@elements),"\n";
	}
}
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
