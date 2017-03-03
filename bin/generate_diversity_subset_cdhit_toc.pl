#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

$| = 1;

# Command line parameters
my ($assemblyDir,$individualMsf,$VERBOSE,$PERCENTRANGE,$MINPERCENT,$assemblyWatered,$individualMsfWatered,$scale,$cdHitFile,$seqcleanFile,$outputDir);
my @options = (
	'a=s',		\$assemblyDir,
	'aw=s',		\$assemblyWatered,
	'i=s',		\$individualMsf,
	'iw=s',		\$individualMsfWatered,
	'v=s',		\$VERBOSE,
	'r=s',		\$PERCENTRANGE,
	'm=s',		\$MINPERCENT,
	'scale=s',	\$scale,
	'c=s',		\$cdHitFile,
	'sc=s',		\$seqcleanFile,
	'o=s',		\$outputDir
);
&GetOptions(@options);
usage() unless 
	defined $assemblyDir
	and defined $individualMsf
	and defined $assemblyWatered
	and defined $individualMsfWatered
	and defined $outputDir
;
$PERCENTRANGE = 1 unless defined $PERCENTRANGE;	# taxon matches within this range will all be considered equal
$MINPERCENT = 80 unless defined $MINPERCENT;	# minimum percent identity required to even consider the taxon label
$scale = 10000 unless defined $scale;

# READING the clustering data to collapse over-split OTUs
# get the length of all of the cleaned up OTU consensus sequences
# assumption: that the isotigs / contigs have been seqcleaned to ensure that primers are removed
$seqcleanFile = join('/',$assemblyDir,'454Isotigs.fna.seqclean') unless defined $seqcleanFile;
die "Error $seqcleanFile does not seem to exist" unless (-e $seqcleanFile);
warn "Reading tig lengths from $seqcleanFile" if $VERBOSE;
my $tigLengths = getLengths($seqcleanFile);


# get the information about which tigs should be substituted based on the cdHit results
# this is so that the clustered (@ 100% identity) results can be reduced to simply the longest sequence
$cdHitFile = $seqcleanFile . '.cd-hit.clstr.report' unless defined $cdHitFile;
die "Error $cdHitFile does not seem to exist" unless (-e $cdHitFile);
warn "Reading tig subsitutions based on $cdHitFile" if $VERBOSE;
my $tigs2substitute = getSubstitutionsFromCDhitClust($cdHitFile,$tigLengths);

# parse the wateredBLAST results
#my $tig2Taxon = parseWateredBLAST($assemblyWatered);	
# parse the wateredBLAST results
my $tig2Taxon;
if (-s $assemblyWatered) {
	warn "Reading wateredBLAST data from $assemblyWatered" if $VERBOSE;
        $tig2Taxon = parseWateredBLAST($assemblyWatered);
}else{
        warn "WARNING $assemblyWatered does not exist so likely everything will be flagged as UNKNOWN";
}


#my $read2Taxon = parseWateredBLAST($individualMsfWatered);
my $toc = generateTocWithSubstitution($assemblyDir,$tigs2substitute);

# do this for all files specified
foreach my $subset (@ARGV){
	my $basename = basename($subset);
	if (!(-e $outputDir)){
		mkdir $outputDir;
	}
	my $output = join('/',$outputDir,$basename);
	die "Error refusing to overwrite $output. This is incase you have multiple files specified on the command line which have the same basename. Check the files in $outputDir and see if you can delete those" if (-e $output);
	open(OUT,">$output") or die "Error opening $output for writting: $!";

	my %diversity;
	my $total = 0;
	
	my %reads2Lookup;
	open(FILE,"<$subset") or die "Error opening $subset: $!";
	while(<FILE>){
		chomp $_;
		my ($read) = split(/\t/,$_);
	#	warn "Reading $read";
		if (defined $toc->{read2Tig}->{$read}){
			$diversity{$toc->{read2Tig}->{$read}}++;
			$total++; # Should ensure that singleton reads are ignored
		}else{
#			$diversity{$read}++;
			$diversity{$read} = 0; # this is a singleton read so blank out its abundance as it should not be used 
			$reads2Lookup{$read}++;
		}
	}
	close FILE or die "Error closing $subset";
	
	#parseSubsetWateredBLAST($individualMsfWatered,\%reads2Lookup,$tig2Taxon);
	if (-s $individualMsfWatered) {
	        parseSubsetWateredBLAST($individualMsfWatered,\%reads2Lookup,$tig2Taxon);
	}else{
	        warn "WARNING $individualMsfWatered does not exist so likely everything will be flagged as UNKNOWN";
	}
	
	
	foreach my $otu (sort{return $diversity{$b} <=> $diversity{$a}} keys %diversity){
		die "Error got $otu somehow when this should have been subsumed" if defined $tigs2substitute->{$otu};
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
	

}
exit 0;
	
sub usage {
	die "Usage: $0 -a assemblyDir -aw assemblyWateredBLAST -i individualMsf -iw individualMsfWateredBLAST -r Range -m minPercent -sc path/to/seqclean/file -c path/to/CD-hit.report -o outputDir ReadFile1 ReadFile2 ReadFile3 ,...";
}

sub getLengths {
        my $file = shift;
        die "Error lost file" unless defined $file;
        my %lengths;
        my $inio = Bio::SeqIO->new(-file => "<$file", -format => 'fasta');
        while(my $seq = $inio->next_seq()){
                $lengths{$seq->display_id()} = $seq->length();
        }
        return \%lengths;
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



sub generateTocWithSubstitution {
	my $dir = shift;
	die "Error lost assembly directory" unless defined $dir;

	my $tigs2substitute = shift; # this is a hashref which takes as a key the name of the OTU consensus which is subsumed by the returning value OTU consensus name

	my %toc;
	my %seenCns;

	my $file = join('/',$dir,'ace.toc');
	warn "Reading assembly location for reads from $file" if $VERBOSE;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	while(<FILE>){
		chomp $_;
		my ($cn,$read) = split(/\s+/,$_);
		warn "Reading data for $cn" if ($VERBOSE and (!defined $seenCns{$cn}));
		$seenCns{$cn}++;
		# if there is a substitution to make based on the clustering results
		# then assign this read to the new isotig
#			$cn = $tigs2substitute->{$cn} if (defined $tigs2substitute->{$cn});
		if (defined $tigs2substitute->{$cn}){
			$cn = $tigs2substitute->{$cn};
		}elsif(!defined $tigLengths->{$cn}){
			warn "$cn data is being skipped - probably got trashed by seqclean";
			next;
		}
		$toc{read2Tig}->{$read} = $cn;
		$toc{tig2Reads}->{$cn}->{$read}++;
	}
	close FILE or die "Error closing $file: $!";

	return \%toc;
}
