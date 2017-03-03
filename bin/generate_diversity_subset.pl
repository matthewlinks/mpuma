#!/usr/bin/perl -w
use strict;
use Getopt::Long;

$| = 1;

# Command line parameters
my ($assemblyDir,$individualMsf,$subset,$VERBOSE,$PERCENTRANGE,$MINPERCENT,$assemblyWatered,$individualMsfWatered,$scale);
my @options = (
	'a=s',		\$assemblyDir,
	'aw=s',		\$assemblyWatered,
	'i=s',		\$individualMsf,
	'iw=s',		\$individualMsfWatered,
	's=s',		\$subset,
	'v=s',		\$VERBOSE,
	'r=s',		\$PERCENTRANGE,
	'm=s',		\$MINPERCENT,
	'scale=s',	\$scale
);
&GetOptions(@options);
usage() unless 
	defined $assemblyDir
	and defined $individualMsf
	and defined $assemblyWatered
	and defined $individualMsfWatered
	and defined $subset
;
$PERCENTRANGE = 1 unless defined $PERCENTRANGE;	# taxon matches within this range will all be considered equal
$MINPERCENT = 80 unless defined $MINPERCENT;	# minimum percent identity required to even consider the taxon label
$scale = 10000 unless defined $scale;

# parse the wateredBLAST results
my $tig2Taxon;
if (-s $assemblyWatered) {
	$tig2Taxon = parseWateredBLAST($assemblyWatered);
}else{
	warn "WARNING $assemblyWatered does not exist so likely everything will be flagged as UNKNOWN";
}
#my $read2Taxon = parseWateredBLAST($individualMsfWatered);
my $toc = generateToc($assemblyDir);

my %diversity;

my $total = 0;

my %reads2Lookup;
open(FILE,"<$subset") or die "Error opening $subset: $!";
while(<FILE>){
	chomp $_;
	my ($read) = split(/\t/,$_);
#	warn "Reading $read";
	$total++;
	if (defined $toc->{read2Tig}->{$read}){
		$diversity{$toc->{read2Tig}->{$read}}++;
	}else{
		$diversity{$read}++;
		$reads2Lookup{$read}++;
	}
}
close FILE or die "Error closing $subset";

if (-s $individualMsfWatered) {
	parseSubsetWateredBLAST($individualMsfWatered,\%reads2Lookup,$tig2Taxon);
}else{
	warn "WARNING $individualMsfWatered does not exist so likely everything will be flagged as UNKNOWN";
}


foreach my $otu (sort{return $diversity{$b} <=> $diversity{$a}} keys %diversity){
	my @matchList;
	if (defined ($tig2Taxon->{$otu})){
		foreach my $match (sort{return $b->{percent} <=> $a->{percent}} @{$tig2Taxon->{$otu}->{matches}}){
			$match->{strand} = '.' unless defined $match->{strand};
			push(@matchList,join(',',$match->{taxon},$match->{percent},$match->{length},$match->{strand}));
		}
	}else{
		push(@matchList,'UNKNOWN');
	}
	my $percent = sprintf("%0.2f",($diversity{$otu} / $total)*100);
	my $scaledLevel = sprintf("%0.2f",(($diversity{$otu} / $total) * $scale));

	print join("\t",$otu,$percent,$scaledLevel,$diversity{$otu},join(';',@matchList)),"\n";
}

exit 0;

sub usage {
	die "Usage: $0 -a assemblyDir -aw assemblyWateredBLAST -i individualMsf -iw individualMsfWateredBLAST -s subsetOfReads -r Range -m minPercent";
}

sub parseWateredBLAST {
	my $wateredBLAST = shift;
	die "Error lost wateredBLAST report" unless defined $wateredBLAST;

	my %query2Taxon;

	open(FILE,"<$wateredBLAST") or die "Error opening $wateredBLAST: $!";
	my $header = <FILE>;
	my $n = 0;
	while(<FILE>){
		chomp $_;
		my ($query,$qDesc,$percentMatch,$lengthMatch,$taxon,$strand) = split(/\t/,$_);
		$n++;
#		warn "$n: $query\n";
		if ($percentMatch < $MINPERCENT){
			# this match is not considered informative
			$taxon = 'UNKNOWN';
		}

		next if $taxon eq 'UNKNOWN'; # skip unknowns so that they don't show up in the SpSn calculations as it would be irrelevant... 
		# they would likely have a Sp of 1 and Sn of 0

		my %match;
		$match{length} = 0;
		$match{percent} = 0;
		$match{percent} = $percentMatch;
		$match{length} = $lengthMatch; 
		$match{taxon} = $taxon;
		$match{strand} = $strand; 

		# figure out whether the new match needs to be added to the list of matches...
		my $addMatch = 0;
		if (!defined @{$query2Taxon{$query}->{matches}}){
			# this is the first match so just add it
			$addMatch = 1;
		}else{
			# see what the current best match is
			my ($bestMatch) = sort{return $b->{percent} <=> $a->{percent}} @{$query2Taxon{$query}->{matches}};
			die "Error lost bestMatch" unless defined $bestMatch;

			# if the new match is bettern than the best match then keep it
			
			if ((($match{percent} - $bestMatch->{percent}) > 0)or(abs($bestMatch->{percent} - $match{percent}) <= $PERCENTRANGE)){
				# it is within $PERCENTRANGE of the best match then definately keep it
				# or
				# it is better than the existing match so keep it for sure
				$addMatch = 1;
			}
		}

		if ($addMatch){
			# add the new match to the list
			push(@{$query2Taxon{$query}->{matches}},\%match) if ($addMatch);

			# need to check and see if the new match causes some of them to be dropped out of the list

			# get a copy of the matches listed by order
			my @sortedMatches =  sort{return $b->{percent} <=> $a->{percent}} @{$query2Taxon{$query}->{matches}};

			my $bestMatch = $sortedMatches[0];

			# delete the previous ones 
			undef @{$query2Taxon{$query}->{matches}};

			# add them back until done or until out of the PERCENTRANGE
			foreach my $m (@sortedMatches){
				if (abs($bestMatch->{percent} - $m->{percent}) <= $PERCENTRANGE){
					push(@{$query2Taxon{$query}->{matches}},$m);
				}else{
					last;
				}
			}
		}
	}
	close FILE or die "Error closing $wateredBLAST: $!";
	return \%query2Taxon;
}

sub parseSubsetWateredBLAST{	
	my $wateredBLAST = shift;
	die "Error lost wateredBLAST report" unless defined $wateredBLAST;
	my $queriesOfNote = shift;
	die "Error lost queries of note" unless defined $queriesOfNote;
	my $query2Taxon = shift;
	die "Error lost query to taxon" unless defined $query2Taxon;

	open(FILE,"<$wateredBLAST") or die "Error opening $wateredBLAST: $!";
	my $header = <FILE>;
	my $n = 0;
	while(<FILE>){
		chomp $_;
		my ($query,$qDesc,$percentMatch,$lengthMatch,$taxon,$strand) = split(/\t/,$_);
		$n++;
		next unless (defined $query2Taxon->{$query});
#		warn "$n: $query\n";
		
		if ($percentMatch < $MINPERCENT){
			# this match is not considered informative
			$taxon = 'UNKNOWN';
		}

		next if $taxon eq 'UNKNOWN'; # skip unknowns so that they don't show up in the SpSn calculations as it would be irrelevant... 
		# they would likely have a Sp of 1 and Sn of 0

		my %match;
		$match{length} = 0;
		$match{percent} = 0;
		$match{percent} = $percentMatch;
		$match{length} = $lengthMatch; 
		$match{taxon} = $taxon;
		$match{strand} = $strand; 

		# figure out whether the new match needs to be added to the list of matches...
		my $addMatch = 0;
		if (!defined @{$query2Taxon->{$query}->{matches}}){
			# this is the first match so just add it
			$addMatch = 1;
		}else{
			# see what the current best match is
			my ($bestMatch) = sort{return $b->{percent} <=> $a->{percent}} @{$query2Taxon->{$query}->{matches}};
			die "Error lost bestMatch" unless defined $bestMatch;

			# if the new match is bettern than the best match then keep it
			
			if ((($match{percent} - $bestMatch->{percent}) > 0)or(abs($bestMatch->{percent} - $match{percent}) <= $PERCENTRANGE)){
				# it is within $PERCENTRANGE of the best match then definately keep it
				# or
				# it is better than the existing match so keep it for sure
				$addMatch = 1;
			}
		}

		if ($addMatch){
			# add the new match to the list
			push(@{$query2Taxon->{$query}->{matches}},\%match) if ($addMatch);

			# need to check and see if the new match causes some of them to be dropped out of the list

			# get a copy of the matches listed by order
			my @sortedMatches =  sort{return $b->{percent} <=> $a->{percent}} @{$query2Taxon->{$query}->{matches}};

			my $bestMatch = $sortedMatches[0];

			# delete the previous ones 
			undef @{$query2Taxon->{$query}->{matches}};

			# add them back until done or until out of the PERCENTRANGE
			foreach my $m (@sortedMatches){
				if (abs($bestMatch->{percent} - $m->{percent}) <= $PERCENTRANGE){
					push(@{$query2Taxon->{$query}->{matches}},$m);
				}else{
					last;
				}
			}
		}
	}
	close FILE or die "Error closing $wateredBLAST: $!";
}


sub generateToc {
	my $dir = shift;
	die "Error lost assembly directory" unless defined $dir;

	my %toc;

	my $aceDir = join('/',$dir,'ace');
	opendir(DIR,$aceDir) or die "Error opening $aceDir for reading: $!";
	my $num = 0;
	while(my $file = readdir(DIR)){
		next if ($file =~ /^\./);
		next unless ($file =~ /^(isotig\d+)\.ace$/);
		my $cn = $1;
		my $aceFile = join('/',$aceDir,$file);
		$num++;
		warn "$num: Reading data from $aceFile" if $VERBOSE;
		open(FILE,"<$aceFile") or die "Error opening $aceFile: $!";
		while(<FILE>){
			chomp $_;
			if ($_ =~ /^RD/){
				my @parts = split(/\s+/,$_);
				my $read = $parts[1];
				while($read =~ /^(.*)\..*$/){
					$read = $1;
				}
				$toc{read2Tig}->{$read} = $cn;
				$toc{tig2Reads}->{$cn}->{$read}++;
			}else{
				if($_ =~ /^RD/){
					die "Error parsing Read line [$_] from $aceFile";
				}else{
					next;
				}
			}
		}
		close(FILE) or die "Error closing $aceFile: $!";
	}
	closedir(DIR) or die "Error closing $aceDir: $!";
	return \%toc;
}
