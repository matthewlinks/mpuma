#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

my ($output,$cdhit_file);

my @options = (
	'o=s',	\$output,
	'c=s',	\$cdhit_file
);
&GetOptions(@options);
usage() unless defined $cdhit_file;


#my $lengths = getLengths($fasta);
my ($clusters,$lengths) = parseCDhitESTClust($cdhit_file);

open(STDOUT,">$output") or die "Error pointing STDOUT to $output: $!" if defined $output;

print join("\t",'Name of longest sequence','Length of longest sequence','# members','CSV of members within this cluster'),"\n";
foreach my $id (sort { return @{$clusters->{$b}} <=> @{$clusters->{$a}}} keys %{$clusters} ){
	my $number = @{$clusters->{$id}};
	print join("\t",$id,$lengths->{$id},$number,join(',',@{$clusters->{$id}})),"\n";
}

exit 0;

sub usage {
	die "Usage: $0 -c cd-hit-est.clstr";
}

sub parseCDhitESTClust {
        my $cdhit = shift;
        die "Error lost blastclust file" unless defined $cdhit;
#        my $tigLength = shift;
#        die "Error lost contig lengths" unless defined $tigLength;

	# Read in all the clusting info
	my (%clusters,%lengths,$representative,$cl,@members);
        open(FILE,"<$cdhit") or die "Error opening $cdhit for reading: $!";
	#check that the first line is a Cluster
	my $line = <FILE>;
	chomp $line;
	if ($line =~ /^\>Cluster\s+(\d+)/){
		$cl = $1;
	}else{
		die "Error expecting Cluster line but got [$line]";
	}
        while(<FILE>){
                chomp $_;
		if ($_ =~ /^\>Cluster\s+(\d+)/){
			my $newCl = $1;
			die "Error I don't have a representative for $cl" unless defined $representative;
			die "Error got multiple clusters defined for $representative" if (defined $clusters{$representative});
			my @tmp_members = @members;
			$clusters{$representative} = \@tmp_members;
			undef @members;
			undef $representative;
			$cl = $newCl;
#		}elsif($_ =~ /\d+\s+(\d+)(nt|aa),\s\>(\w+)\.\.\.\s\*/){
#			die "Error got multiple representatives defined for $cl" if defined $representative;
#			$representative = $3;
#			die "Error multiple lengths reported for $representative" if defined $lengths{$representative};
#			$lengths{$representative} = $1;
#		}elsif($_ =~ /\d+\s+(\d+)(nt|aa),\s\>(\w+)\.\.\.\sat/){
#			die "Error multiple lengths reported for $3" if defined $lengths{$3};
#			$lengths{$3} = $1;
#			push(@members,$3);
# if we get Sanger data then it looks like
#0	727nt, >SH0807C-T7-001_C04_14AUG2012_028.ab1... *

		}elsif($_ =~ /\d+\s+(\d+)(nt|aa),\s\>([\w\_\-\.]+)\.\.\.\s\*$/){
			die "Error got multiple representatives defined for $cl" if defined $representative;
			$representative = $3;
			die "Error multiple lengths reported for $representative" if defined $lengths{$representative};
			$lengths{$representative} = $1;
		}elsif($_ =~ /\d+\s+(\d+)(nt|aa),\s\>(\w+)\.\.\.\sat/){
			die "Error multiple lengths reported for $3" if defined $lengths{$3};
			$lengths{$3} = $1;
			push(@members,$3);
		}else{
			die "Error cannot parse [$_]";
		}
	}
	die "Error I don't have a representative for $cl" unless defined $representative;
	die "Error got multiple clusters defined for $representative" if (defined $clusters{$representative});
	my @tmp_members = @members;
	$clusters{$representative} = \@tmp_members;
#	$clusters{$representative} = @members;
        close FILE or die "Error closing $cdhit: $!";
        return (\%clusters,\%lengths);
}

