#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

$| = 1;

# Command line parameters
my ($VERBOSE,$PERCENTRANGE,$MINPERCENT,$watered,$scale,$cdHitFile,$outputFile,$inputFile,$tocFile);
my @options = (
	'w=s',		\$watered,
	'v=s',		\$VERBOSE,
	'r=s',		\$PERCENTRANGE,
	'm=s',		\$MINPERCENT,
	'scale=s',	\$scale,
	'c=s',		\$cdHitFile,
	'i=s',		\$inputFile,
	'o=s',		\$outputFile,
	't=s',		\$tocFile
);
&GetOptions(@options);
usage() unless defined $tocFile
	and defined $inputFile
	and defined $outputFile
;
$PERCENTRANGE = 1 unless defined $PERCENTRANGE;	# taxon matches within this range will all be considered equal
$MINPERCENT = 80 unless defined $MINPERCENT;	# minimum percent identity required to even consider the taxon label
$scale = 10000 unless defined $scale;

# make sure we aren't overwritting data
die "Error refusing to overwrite $outputFile." if (-s $outputFile);

# figure out what substitutions are to be made if there is a CD-hit clustering report provided
# this is usually provided when looking for clusters of synonymous DNA changes (i.e. non-redundant OTU in the protein space)
my $OTUs2sub;
$OTUs2sub = getSubstitutionsFromCDhitClust($cdHitFile) if (defined $cdHitFile);

my ($readToNrOTU,$nrOTUtoRead) = parse_toc($tocFile,$OTUs2sub);

# parse the wateredBLAST results
my $tig2Taxon;
if (-s $watered) {
	warn "Reading wateredBLAST data from $watered" if $VERBOSE;
        $tig2Taxon = parseWateredBLAST($watered);
}else{
        warn "WARNING $watered does not exist so likely everything will be flagged as UNKNOWN";
}


# calculate the actual diversity
my %diversity;
my $total = 0;

my %reads2Lookup;
open(FILE,"<$inputFile") or die "Error opening $inputFile: $!";
while(<FILE>){
	chomp $_;
	my ($read) = split(/\t/,$_);
	my $nrOTU = $readToNrOTU->{$read};
	if ($nrOTU eq 'Singleton'){
		# sink all the Singleton data so that they don't appear anywhere
#		$diversity{$read} = 0; # Singletons don't count
#		$reads2Lookup{$read}++;
	}elsif($nrOTU eq 'None'){
		# trashed by seqclean
		# these reads don't count
	}else{
		$diversity{$nrOTU}++;
		$total++;		# reads which are NOT singletons count for total abundance
	}
}
close FILE or die "Error closing $inputFile";

# Generate the output
open(OUT,">$outputFile") or die "Error opening $outputFile for writting: $!";
foreach my $otu (sort{return $diversity{$b} <=> $diversity{$a}} keys %diversity){
	die "Error got $otu somehow when this should have been subsumed" if defined $OTUs2sub->{$otu};
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

	print OUT join("\t",$otu,$percent,$scaledLevel,$diversity{$otu},join(';',@matchList)),"\n";
}


exit 0;
	
sub usage {
	die "Usage: $0 -t total.toc -i inputReadFile -o outputDiversityFile";
}

sub parseWateredBLAST {
	my $wateredBLAST = shift;
	die "Error lost wateredBLAST report" unless defined $wateredBLAST;

	my %query2Taxon;

	open(FILE,"<$wateredBLAST") or die "Error opening $wateredBLAST: $!";
#	my $header = <FILE>;
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

sub getSubstitutionsFromCDhitClust {
        my $cdHit = shift;
        die "Error lost CD-hit file" unless defined $cdHit;

        my (%IDs2skip);
        open(FILE,"<$cdHit") or die "Error opening $cdHit for reading: $!";
	my $header = <FILE>;
        while(<FILE>){
                chomp $_;
		my ($largest_id,$length,$num_members,$csv_members) = split(/\t/,$_);
		my @smaller_ids = split(/\,/,$csv_members);
                foreach my $id (@smaller_ids){
                        $IDs2skip{$id} = $largest_id;
                }
        }
        close FILE or die "Error closing $cdHit: $!";
        return \%IDs2skip; # a hashref of which substitutions should be made
}

sub parse_toc {
	my $file = shift;
	die "Error lost toc file" unless defined $file;
	my $tigs2skip = shift;

	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $header = <FILE>;
	my %read2OTU;
	my %OTU2read;
	while(<FILE>){
		chomp $_;
		my ($read,$nrOTU,$rationale) = split(/\,/,$_);
		# make the substitution if it is defined
		$nrOTU = $OTUs2sub->{$nrOTU} if ((defined $OTUs2sub) and (defined $OTUs2sub->{$nrOTU}));
		$read2OTU{$read} = $nrOTU;
		push(@{$OTU2read{$nrOTU}},$read);
	}
	close FILE or die "Error closing $file: $!";
	return (\%read2OTU,\%OTU2read);
}


