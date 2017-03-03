#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use IPC::Open2;
use File::Basename;
use IO::Compress::Gzip;
use IO::Uncompress::Gunzip;
use Bio::SeqIO;

my ($hostGenome,$cpu,$outDir,$compressed,$fastqType);
my @options = (
	'h=s',	\$hostGenome,
	'c=s',	\$cpu,
	'o=s',	\$outDir,
	'z=s',	\$compressed,
	't=s',	\$fastqType
);
&GetOptions(@options);
$compressed = 0 unless defined $compressed;
$fastqType = 'illumina' unless defined $fastqType;
my %supportedFastqTypes = (
	'sanger'	=> 1,
	'solexa'	=> 1,
	'illumina'	=> 1
);
die "Error $fastqType is not a supported type of fastq encoding"
unless defined $supportedFastqTypes{$fastqType};

my $seqIOtype = join('-','fastq',$fastqType);

die "Usage: $0 -h host_genome -c cpuNum -o outDir -t illumina|solexa|sanger -z 0 file0 file1 file2 file3 ..." 
unless(
	defined $hostGenome
	and defined $cpu
	and defined $outDir
);

foreach my $file (@ARGV){
	my $basename = basename($file);
	my $output = join('/',$outDir,$basename);
	die "Error $output already exists...refusing to overwrite" if (-e $output);
	my $noHitIDs = get_ids_from_bt2($file,$hostGenome,4,$cpu);
	my ($in,$out,$inio,$outio);	
	if ($compressed){
		$in = IO::Uncompress::Gunzip->new($file) or die "Error opening $file for reading: $!";
		$inio = Bio::SeqIO->new(-Fh => $in, -format => $seqIOtype);
		$out = IO::Compress::Gzip->new($output) or die "Error opening $output for reading: $!";
		$outio = Bio::SeqIO->new(-Fh => $out, -format => $seqIOtype);
	}else{
		$inio = Bio::SeqIO->new(-file => "<$file", -format => $seqIOtype);
		$outio = Bio::SeqIO->new(-file => ">$output", -format => $seqIOtype);
	}

	# write only those reads which fail to map to the host genome
	while(my $seq = $inio->next_seq()){
		next unless defined $noHitIDs->{$seq->display_id};
		$outio->write_seq($seq);
	}

	# undef all the file stuff just in case something stays open
	undef $inio;
	undef $in;
	undef $outio;
	undef $out;
}

exit 0;

sub get_ids_from_bt2 {
	my $query = shift;
	die "Error lost query" unless defined $query;
	my $base = shift;
	die "Error lost base" unless defined $base;
	my $flag = shift;
	die "Error lost flag" unless defined $flag;
	my $cpu = shift;
	$cpu = 2 unless defined $cpu;

	my %ids;
	my $pid = open2(\*RDR,\*WTR,'bowtie2',
		'--local',
		'-t',
		'-p',	$cpu,
		$base,
		'-U',	$query,
		'--no-head') or die "Error calling bowtie2: $!";
	close WTR or die "Error closing writter to bowtie2 process: $!";
	while(<RDR>){
		chomp $_;
		my ($id,$f) = split(/\t/,$_);
		$ids{$id}++ if ($f == $flag); 
	}
	close RDR or die "Error closing reader from bowtie2 process: $!";
	wait;
	return \%ids;
}
