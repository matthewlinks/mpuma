#!/usr/bin/perl -w
use strict;
use Getopt::Long;

$| = 1;

# Command line parameters
my ($assemblyDir,$individualMsf,$skip,$v,$database,$VERBOSE,$contigFileName,$PERCENTRANGE,$MINPERCENT,$output,$individualMsfWatered,$contigsWatered);
my @options = (
	'a=s',		\$assemblyDir,
	'i=s',		\$individualMsf,
	'iw=s',		\$individualMsfWatered,
	's=s',		\$skip,
	'v=s',		\$v,
	'd=s',		\$database,
	'verbose=s',	\$VERBOSE,
	'c=s',		\$contigFileName,
	'cw=s',		\$contigsWatered,
	'r=s',		\$PERCENTRANGE,
	'm=s',		\$MINPERCENT,
	'o=s',		\$output,
);
&GetOptions(@options);
usage() unless 
	defined $assemblyDir
	and defined $individualMsf
	and defined $database
;
$skip = 0 unless defined $skip;
$PERCENTRANGE = 1 unless defined $PERCENTRANGE;	# taxon matches within this range will all be considered equal
$MINPERCENT = 80 unless defined $MINPERCENT;	# minimum percent identity required to even consider the taxon label

if (defined $output){
	open(STDOUT,">$output") or die "Error opening STDOUT to point at $output: $!";
}

$contigFileName = '454Isotigs.fna' unless defined $contigFileName;
my $contigs = join('/',$assemblyDir,$contigFileName);
# skip the wateredBLAST step if instructed to do so
if ($skip){
	$contigsWatered = $contigs . '.wateredBLASTN'  unless defined $contigsWatered;
	$individualMsfWatered = $individualMsf . '.wateredBLASTN' unless defined $individualMsfWatered;
}else{
	$contigsWatered = waterBLAST($contigs,$database,$v);
	$individualMsfWatered = waterBLAST($individualMsf,$database,$v);
}

# parse the wateredBLAST results
my $tig2Taxon = parseWateredBLAST($contigsWatered);	
my $read2Taxon = parseWateredBLAST($individualMsfWatered);
my $toc = generateToc($assemblyDir);

# calculate the Spesificity and Sensitivity metrics for each tig
print join("\t",'Tig','Specifity','Sensitivity','Residual Distance','Number of Reads','Matches'),"\n";

#foreach my $tig (keys %{$tig2Taxon}){
foreach my $tig (sort {
		return $a cmp $b
		}keys %{$toc->{tig2Reads}}){
	
	warn "TIG $tig";
	my ($Sp,$Sn) = ('N/A','N/A');
	my $residual = 'N/A';
	my @matchList;
	if (defined $tig2Taxon->{$tig}){
		# this is not and UNKNOWN so we can try to assess the Sp and Sn for it
	
		warn "Calculate SpSn for $tig";
		($Sp,$Sn) = calculateSpSn($tig,$tig2Taxon,$read2Taxon,$toc);
		$residual = calculateResidual($Sp,$Sn); # Calculate the residual distance to perfection 1,1
		warn "$tig: (Sp,Sn) = ($Sp,$Sn)\n";

		foreach my $match (sort{return $b->{percent} <=> $a->{percent}} @{$tig2Taxon->{$tig}->{matches}}){
			$match->{strand} = '.' unless defined $match->{strand};
			push(@matchList,join(',',$match->{taxon},$match->{percent},$match->{length},$match->{strand}));
		}
	}

	my $numReads = keys %{$toc->{tig2Reads}->{$tig}};
	if (@matchList > 0){
		print join("\t",$tig,$Sp,$Sn,$residual,$numReads,join(";",@matchList)),"\n";
	}else{
		print join("\t",$tig,$Sp,$Sn,$residual,$numReads,'UNKNOWN'),"\n";
	}
}


#foreach my $tig (keys %{$toc->{tig2Reads}}){
#	print "$tig:";
#	my @reads = keys %{$toc->{tig2Reads}->{$tig}};
#	my $num = @reads;
##	print "\t",join(',',@reads),"\n";
#	print "\t",$num,"\n";
#
#}

exit 0;

sub calculateResidual {
	my $Sp = shift;
	die "Error lost Specificity" unless defined $Sp;
	my $Sn = shift;
	die "Error lost Sensitivity" unless defined $Sn;

	if ($Sp ne 'N/A'){
		die "Error Specificity is > 1" unless ($Sp <= 1);
		die "Error Specificity is < 0" unless ($Sp >= 0);
	}
	if ($Sn ne 'N/A'){
		die "Error Sensitivity is > 1" unless ($Sn <= 1);
		die "Error Sensitivity is < 0" unless ($Sn >= 0);
	}

	my $residual = 0;

	if ($Sp == 1){
		$residual = sprintf("%0.2f",1-$Sn);
	}elsif($Sn == 1){
		$residual = sprintf("%0.2f",1-$Sp);
	}else{
		# solve for hypoteneuse
		my $aSquared = (1-$Sp)*(1-$Sp);
		my $bSquared = (1-$Sn)*(1-$Sn);
		$residual = sprintf("%0.2f",sqrt($aSquared + $bSquared));
	}
	return $residual;
}

sub calculateSpSn {
	my $tig = shift;
	die "Error lost tig" unless defined $tig;
	my $tig2Taxon = shift;
	die "Error lost tig2Taxon" unless defined $tig2Taxon;
	my $read2Taxon = shift;
	die "Error lost read2Taxon" unless defined $read2Taxon;
	my $toc = shift;
	die "Error lost Table of Contents" unless defined $toc;
	
	# for simplicity create a hash containing the taxons which are considered as positives for this contig
	my %tigTaxon;
	foreach my $match (@{$tig2Taxon->{$tig}->{matches}}){
		$tigTaxon{$match->{taxon}}++;
	}
	
	my ($TN,$TP,$FN,$FP) = (0,0,0,0);

	# calculate True Positives (TP)
	# this is the number of reads in this tig which agree with the classification of the tig itself
	# calculate False Positives (FP)
	# this is the number of reads in this tig which do NOT agree with the classification of the tig itself
	foreach my $read (keys %{$toc->{tig2Reads}->{$tig}}){
		# does the read agree with the tig
		my $agree = 0;
		foreach my $m (@{$read2Taxon->{$read}->{matches}}){
			if ($tigTaxon{$m->{taxon}}){
				$agree = 1;
				last;
			}
		}
		if ($agree){
			# if it agrees with the Taxon for this tig then this is a true positive (TP)
			$TP++;
		}else{
			# it disagrees with the taxon for this tig then this read is a false positive (FP)
			$FP++;
		}
	}
		
	# calculate the True negatives (TN)
	# this is the number of reads in other tigs (with other taxon matches)
	# calculate the False negatives (FN) this is the number of reads in other tigs which should have been in this tig
	foreach my $otherTig (keys %{$tig2Taxon}){
		next if $otherTig eq $tig; # skip self
		my $tigAgree = 0;
		foreach my $m (@{$tig2Taxon->{$otherTig}->{matches}}){
			if ($tigTaxon{$m->{taxon}}){
				$tigAgree = 1;
				last;
			}
		}
		if ($tigAgree){
			# then these two tigs are matching similar or the same taxon could be 5' vs. 3' end of the target or highly related 
			# (refinement of which would take a look at the % ID relative to the taxon etc.)
			# either way this otherTig cannot be used to calculate TN or FN
			next;
		}else{
			# look through the reads for this otherTig for reads which look like $tigTaxon 
			foreach my $read (keys %{$toc->{tig2Reads}->{$otherTig}}){
				# does the read agree with $tigTaxon
		                my $agree = 0;
        		        foreach my $readM (@{$read2Taxon->{$read}->{matches}}){
  		                      if ($tigTaxon{$readM->{taxon}}){
						$agree = 1;
						last;
		                        }
		                }
		                if ($agree){
		                        # if it agrees with the tigTaxon for this tig then this is a false negative (FN)
		                        $FN++;
		                }else{
                		        # it disagrees with the taxon for this tig then this read is a true negative (FP)
                		        $TN++;
               		 	}
			}
		}
        }

	warn "TN,TP,FP,FN = $TN,$TP,$FP,$FN";

	my $Sp = 'N/A';
	$Sp = sprintf("%0.2f",($TN/($TN+$FP))) if (($TN + $FP) != 0);

	my $Sn = 'N/A';
	$Sn = sprintf("%0.2f",($TP/($TP+$FN))) if (($TP+$FN) != 0);
	return ($Sp,$Sn);
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
		warn "$n: $query\n";
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
sub usage {
	die "Usage: $0 -a assemblyDir -i individualMsf -s skipAndUseThisMsf -v 5 -d database -o output";
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

sub waterBLAST {
	my $input = shift;
	die "Error lost input file" unless defined $input;
	my $database = shift;
	die "Error lost database" unless defined $database;
	my $v = shift;
	$v = 5 unless defined $v;
	my $output = $input . '.wateredBLASTN';
	
	system("$ENV{METACPN}/bin/watered_blast.pl",
		'-i',	$input,
		'-d',	$database,
		'-o',	$output,
		'-v',	$v,
		'-p',	'blastn'
	) == 0 or die "Error calling waterBLAST for $input: $?";
	return $output;
}

