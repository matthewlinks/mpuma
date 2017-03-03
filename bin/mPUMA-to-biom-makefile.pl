#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use JSON;
use Bio::SeqIO;
use IPC::Open2;
use POSIX qw(strftime);

my (
$OTUfastaFile,$OTUwbFile,$geneTarget,
$totalToc,$chimeraBoolean,$assemblyMethod,$mappingMethod,
$isProtein,$outputFile,$libraryMetadataFile,$libraryMetadataDelim,
$VERBOSE,$prettyPrint,$OTUtree,$OTUproteinWbFile,$librarySpec,$clusteringReport,
$forceFasta,
);

my $VERSION = '1.0';
my @options = (
	'a=s',	\$assemblyMethod,
	'c=s',	\$chimeraBoolean,
	'cr=s',	\$clusteringReport,
	'm=s',	\$mappingMethod,
	'of=s',	\$OTUfastaFile,
	'ff=s',	\$forceFasta,
	'ow=s',	\$OTUwbFile,
	'pw=s',	\$OTUproteinWbFile,
	'ot=s',	\$OTUtree,
	'g=s',	\$geneTarget,
	'p=s',	\$isProtein,
	'ls=s',	\$librarySpec,
	'lm=s',	\$libraryMetadataFile,
	'ld=s',	\$libraryMetadataDelim,
	'o=s',	\$outputFile,
	'v=s',	\$VERBOSE,
	'pp=s',	\$prettyPrint,
	'toc=s',	\$totalToc,
);
&GetOptions(@options);

$VERBOSE			= 0 unless defined $VERBOSE;
$prettyPrint			= 0 unless defined $prettyPrint;
$libraryMetadataDelim		= "\t" unless defined $libraryMetadataDelim;
$forceFasta			= 0 unless defined $forceFasta;

die <<EOF
Usage: $0 -a gsAssembler -c 0 -m bowtie2 -g cpn60 -ls /path/to/library.spec -toc /path/to/total.toc

Required: 
	-a  gsAssembler			what assembly approach did mPUMA use
	-m  bowtie2			what mapping technique did mPUMA use
	-c  0				boolean for whether mPUMA's C3 chimera check was used as a basis to EXCLUDE OTUs
	-g  cpn60			this is a simple label for what this experiment was trying to use as a target barcode
	-ls library.spec		the library to read file specification
	-toc total.toc			the total table of contents file. Built by mPUMA to record where each read has been assigned

Options / override-able:
	-of /path/to/OTU.fasta		this is the location of the OTU sequences to pull their sequences and include them from 
	-ff 0|1				this is a boolean used to enforce that wether ONLY the OTUs in the FASTA file are reported
	-ow /path/to/OTU.wateredBLAST	this is the location of a wateredBLAST output file where the information from the best hit will be pulled
	-ot /path/to/OTU.newick		text file which contains a treefile you would like to include
	-p  0				a boolean on whether the OTUs are based on protein sequences
	-lm /path/to/lib/metadata	additional metadata (e.g. experimental factors) which can be included
	-ld '\\t'			the delimiter used in the metadata file
	-o  /path/to/output.biom	where to write the output (default is STDOUT)
	-pp 0				boolean for whether to pretty print the JSON
EOF
unless(
	defined $assemblyMethod 
	and defined $chimeraBoolean
	and defined $mappingMethod
	and defined $geneTarget
	and defined $librarySpec
	and defined $totalToc
);
# start building the object
my %obj;

# setup the header information
warn "Setting up header info" if $VERBOSE;
set_header_info(\%obj);

# add the assembly method
warn "Setting the assembly method" if $VERBOSE;
$obj{'mPUMA'}->{'assemblyMethod'} = $assemblyMethod;

# add the mapping method
warn "Setting mapping method" if $VERBOSE;
$obj{'mPUMA'}->{'mappingMethod'} = $mappingMethod;

# add the chimera boolean
warn "Setting the chimera usage" if $VERBOSE;
$obj{'mPUMA'}->{'C3chimeraUsage'} = $chimeraBoolean;

# add the gene target
warn "Setting the gene target" if $VERBOSE;
$obj{'mPUMA'}->{'geneTarget'} = $geneTarget;

######################### Library section

# need a generic hash to hold metadata and analysis info on the libraries
my %genericLibs;
warn "Getting all the files from $librarySpec" if $VERBOSE;
my $lib2files = get_files_by_library($librarySpec);
warn "Getting IDs for each library" if $VERBOSE;
my $lib2ids = get_ids_by_library($lib2files);

######################### Add Generic metadata utility added here....
warn "Checking for additional outputs" if $VERBOSE;
parse_additional_library_metadata(\%genericLibs,$libraryMetadataFile,$libraryMetadataDelim) if defined $libraryMetadataFile;

# add any other info about or derived from libraries here

######################### OTU section

# figure out which reads are in which otus and allow going from otu => reads or read => otus
# notice the chance for bidirectional 1 to many relationships
# this is just to support cases where a read could have equivalently good matches
my ($id2otu,$otu2id) = get_id_otu_map($totalToc);

### should have a condensation possible here if there is a clustering report specified 
condense_based_on_clustering($clusteringReport,$id2otu,$otu2id) if defined ($clusteringReport);

# need a generic hash to hold metadata and analysis info on the OTU
my %genericOTUs;

# read OTU sequences
warn "Checking for OTU sequence data" if $VERBOSE;
parse_OTU_sequences(\%genericOTUs,$OTUfastaFile) if (defined $OTUfastaFile);

# read wb results
warn "Checking for wateredBLAST results" if $VERBOSE;
parse_OTU_wateredBLAST(\%genericOTUs,$OTUwbFile) if (defined $OTUwbFile);
parse_OTU_wateredBLAST(\%genericOTUs,$OTUproteinWbFile,'Protein') if (defined $OTUproteinWbFile);

warn "Checking for OTU tree to include" if $VERBOSE;
$obj{'OTUtree'} = read_tree($OTUtree) if (defined $OTUtree);

######################## Diversity section

# this is where the actual OTU abundance info is sourced

# If we have been instructed to enforce that the ONLY OTUs reported are from the FASTA 
# N.B. if this is happening then it is likley there will be a disconnect between the totalReadCount|subsampledReadCount and the # of obs 
# reported in the abundance table
my $OTUsInFasta;
if ($forceFasta){
	warn "Enforcing agreement with FASTA [$forceFasta]";
	$OTUsInFasta = getIDs($OTUfastaFile);
	# here we screen the list of OTUs and undef any that were not in the FASTA file
	# If the Fasta was not specified then this should result in no data table
	foreach my $otu (keys %genericOTUs){
		undef $genericOTUs{$otu} unless defined $OTUsInFasta->{$otu};
	}

}

# library stuff to pull read counts and diversity data 
warn "Finding libraries and setting reads counts" if $VERBOSE;
my %abundanceData;

# get all the ids / library
foreach my $lib (keys %{$lib2ids}){
	my %otuAb;
	warn "Calculating abundance data for $lib" if $VERBOSE;
	# how many reads are there
	my $numReads = @{$lib2ids->{$lib}};
	$genericLibs{$lib}->{'metadata'}->{'totalReadCount'}		= $numReads;
	
	# in this case there is no subsampling so set that up
	$genericLibs{$lib}->{'metadata'}->{'subsampledReadCount'}	= $numReads;

	# recored the per otu abundance for this library
	foreach my $read (@{$lib2ids->{$lib}}){
		foreach my $otu (keys %{$id2otu->{$read}}){
			if ($forceFasta){
				next unless defined $OTUsInFasta->{$otu}; # If we are forcing the FASTA to be the limit on the OTUs then we need to skip counting the missing onts
			}
			$otuAb{$otu} += $id2otu->{$read}->{$otu}; # this should be +1 for each read but in case we want to overload it later pass the value through here
		}
	}

	# for the library place the abundance into the overall abundance table
	foreach my $otu (keys %otuAb){
		$abundanceData{$otu}->{$lib} = ($otuAb{$otu} + 0); #### WARNING the + 0 or *1 here is important to convince the PERL interpreter that this is a numeric type and avoid putting quotes around it when outputting in JSON
                $genericOTUs{$otu}->{'analysis'}->{'prevalence'}++;     # count its occurrance based on having been seen in this $lib. N.B. this has no assurance that the abundance was > 0 just that there was a line for it...	
	}

}

#########################################################################
# All metadata and analysis for Libs and OTUs MUST precede this section #
#########################################################################

# Determine an Order for the OTUs and the Libraries
my @sortedOTU = sort {
	return $genericOTUs{$b}->{'analysis'}->{'prevalence'} <=> $genericOTUs{$a}->{'analysis'}->{'prevalence'}
} keys %genericOTUs;	# sorted by prevalence

my @sortedLibs = sort {
	return $a cmp $b
} keys %genericLibs; # sorted libs on their name


# figure out what the dimensions of the OTU table will be
my $numberOTU = @sortedOTU;
my $numberLibs = @sortedLibs;

# set the dimensionality of the shape key
warn "Setting the dimensionality" if $VERBOSE;
$obj{shape} = [$numberOTU,$numberLibs];

# add the rows and the metadata for the OTU
warn "Setting the OTU information" if $VERBOSE;
foreach my $OTU (@sortedOTU){
	my $O = $genericOTUs{$OTU};
	$O->{id} = $OTU;
	push(@{$obj{rows}},$O);
}

# add the columns
warn "Setting the Library information" if $VERBOSE;
foreach my $lib (@sortedLibs){
	my $L = $genericLibs{$lib};
	$L->{id} = $lib;
	push(@{$obj{columns}},$L);
}

# add the data itself
warn "Setting the abundance table" if $VERBOSE;
foreach my $OTU (@sortedOTU){
	my @vals;
	foreach my $lib (@sortedLibs){
		my $v = 0; 
		$v = $abundanceData{$OTU}->{$lib} if (defined $abundanceData{$OTU}->{$lib});
		push(@vals,$v);
	}
	push(@{$obj{data}},\@vals);
}

# print to a file if that is what was called for
open(STDOUT,">$outputFile") or die "Error opening $outputFile for writting: $!" 
	if defined $outputFile;

# convert to JSON
my $jsonOutput;

if ($prettyPrint){
	$jsonOutput = to_json(\%obj, { ascii => 1, pretty => 1 } );
}else{
	$jsonOutput = to_json(\%obj);
}

# print the object
print $jsonOutput;

### done
exit 0;

sub set_header_info {
	my $obj = shift;
	if (!defined $obj){
		my %new;
		$obj = \%new;
	}
	
	# simple header stuff...
	$obj->{'id'}			= '';
	$obj->{'format'}		= 'Biological Observation Matrix 0.9.1-dev';
	$obj->{'format_url'}		= 'http://biom-format.org/documentation/format_versions/biom-1.0.html';
	$obj->{'type'}			= 'OTU table';
	$obj->{'generated_by'}		= 'mPUMA';
	$obj->{'mPUMA'}->{'version'}	= $VERSION;
	$obj->{'date'}			= get_date();
	$obj->{'matrix_type'}		= 'dense';
	$obj->{'matrix_element_type'}	= 'float';
}

sub get_date{
        my $now = time();
        return strftime("%Y-%m-%dT%H:%M:%S", localtime($now));
}

sub parse_OTU_sequences {
	my $obj = shift;
	die "Error obj" unless defined $obj;
	my $file = shift;
	die "Error lost file" unless defined $file;
	my %OTU2Seq;	
	my $inio = Bio::SeqIO->new(-file => "<$file", -format => 'fasta');
	while(my $seq = $inio->next_seq()){
		$obj->{$seq->display_id()}->{'metadata'}->{'sequence'} = $seq->seq();
	}
}

sub parse_OTU_wateredBLAST {
	my $obj = shift;
	die "Error lost obj" unless defined $obj;
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $string = shift;
	$string = 'Nucleotide' unless defined $string;

	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my %results;
	while(<FILE>){
		chomp $_;
		my (%hit,$id);
		($id,$hit{desc},$hit{percent},$hit{'length'},$hit{match},$hit{orientation}) = split(/\t/,$_);
		next unless ($hit{percent} >= 55); # winnow the garbage out this is based on observations of cpn60 data
		next if ((defined $results{$id})and($results{$id}->{percent} >= $hit{percent}));
		$results{$id} = \%hit;
	}
	close FILE or die "Error closing $file: $!";

	foreach my $id (keys %results){
		$obj->{$id}->{'analysis'}->{'wateredBLAST'}->{'nearest'.$string.'Match'}		= $results{$id}->{match};
		$obj->{$id}->{'analysis'}->{'wateredBLAST'}->{'nearest'.$string.'Identity'}		= $results{$id}->{percent};
		$obj->{$id}->{'analysis'}->{'wateredBLAST'}->{'nearest'.$string.'MatchLength'}		= $results{$id}->{'length'};
		$obj->{$id}->{'analysis'}->{'wateredBLAST'}->{'nearest'.$string.'MatchOrientation'}	= $results{$id}->{orientation};
	
	}
}


sub count_lines_in_text_file{
	my $file = shift;
	return 0 unless defined $file;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $count = 0;
	while(<FILE>){
		$count++;
	}
	close FILE or die "Error closing $file: $!";
	return $count;
}


sub parse_additional_library_metadata {
	my $obj = shift;
	die "Error lost file" unless defined $obj;
        my $file = shift;
        die "Error lost file" unless defined $file;
        my $delim = shift;
        die "Error lost delimiter" unless defined $delim;

        open(FILE,'<',$file) or die "Error opening $file for reading: $!";
        my $header = <FILE>;
        chomp $header;
        my @fields = split(/$delim/,$header);
        my $idCol = shift @fields;

        while(<FILE>){
                chomp $_;
                my @values = split(/$delim/,$_);
                my $id = shift @values;
                for(my $i = 0; $i < @values;$i++){
                        die "Error no key for field $i @ line [$_]" unless defined $fields[$i];
			if ($fields[$i] =~ /^\*(.*)$/){
				# pick up on asterix as a notation that this is an experimental factor
				$obj->{$id}->{'metadata'}->{'experimentalFactors'}->{$1} = $values[$i];
			}else{
				$obj->{$id}->{'metadata'}->{$fields[$i]} = $values[$i];
			}
                }
        }
        close FILE or die "Error closing $file: $!";
}


sub read_tree {
	my $file = shift;
	die "Error lost file" unless defined $file;

	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $tree;
	while(<FILE>){
		chomp $_;
		$tree .= $_;
	}
	close FILE or die "Error closing $file: $!";
	return $tree;
}


sub get_ids_by_library {
	my $l2f = shift;
	die "Error lost lib => file map" unless defined $l2f;
	my %lib2ids;
	foreach my $lib (keys (%{$l2f})){
		my @files;
		foreach my $f (keys %{$l2f->{$lib}}){
			push(@files,$f);
		}
		$lib2ids{$lib} = get_ids_from_files(@files);
	}
	return \%lib2ids;
}

sub get_ids_from_files {
	my @ids;
        my $pid = open2(\*IDS_OUT,\*IDS_IN,"$ENV{mPUMA}/bin/dump_ids_from_sff_or_fastq.pl",@_)
                or die "Error calling dump_ids_from_sff_or_fastq.pl ".join(" ",@_).": $!";
        close(IDS_IN) or die "Error closing STDIN to ID stream: $!";
        while(<IDS_OUT>){
		chomp $_;
		my ($id) = split(/\s+/,$_);
		push(@ids,$id);
	}
        wait;
	return \@ids
}

sub get_files_by_library {

	my $spec = shift;
	die "Error lost spec file" unless defined $spec;
	my %lib2file;
	open(SPEC,"<$spec") or die "Error opening $spec for reading: $!";
	while(<SPEC>){
		chomp $_;
		my ($lib,$file) = split(/\t/,$_);
		$lib2file{$lib}->{$file}++;
	}
	close SPEC or die "Error closing $spec: $!";
	return \%lib2file;
}


sub get_id_otu_map {
	my $file = shift;
	die "Error lost file (total.toc)" unless defined $file;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my $header = <FILE>;
	my (%id2otu,%otu2id);
	while(<FILE>){
		chomp $_;
		my ($id,$otu,$rationale) = split(/\,/,$_);
		die "Error lost id" unless defined $id;
		die "Error lost otu" unless defined $otu;
		$id2otu{$id}->{$otu}++;
		$otu2id{$otu}->{$id}++;
	}
	close FILE or die "Error closing $file: $!";
	return (\%id2otu,\%otu2id);
}	

sub condense_based_on_clustering{
	my $reportFile = shift;
	die "Error lost clustering report" unless defined $reportFile;
	my $id2otu = shift;
	die "Error lost id2otu" unless defined $id2otu;
	my $otu2id = shift;
	die "Error lost otu2id" unless defined $otu2id;

	my $otuSubs = get_subs_from_clustering_report($reportFile);

	# make substitions in the otu2id
	foreach my $otu (keys %{$otu2id}){
		next unless defined $otuSubs->{$otu};
		my $subs = $otuSubs->{$otu};
		foreach my $read ($otu2id->{$otu}){
			# add the read to the substitution
			$otu2id->{$subs}->{$read}++;
		}
		# oblate the old otu
		undef $otu2id->{$otu};
	}

	# make the substitutions in id2otu
	foreach my $read (keys %{$id2otu}){
		foreach my $otu (keys %{$id2otu->{$read}}){
			if (defined $otuSubs->{$otu}){
				# note that the read is mapping to the substitution
				$id2otu->{$read}->{$otuSubs->{$otu}}++;
				# remove the count to the otu being substituted away
				$id2otu->{$read}->{$otu}--;
				# clean up the ref if its 0 
				undef $id2otu->{$read}->{$otu} if ($id2otu->{$read}->{$otu} == 0);
			}
		}
	}

}

sub get_subs_from_clustering_report {
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

sub getIDs {
	my $file = shift;
	die "Error lost file" unless defined $file;
	open(FILE,"<$file") or die "Error opening $file for reading: $!";
	my %IDs;
	while(<FILE>){
		chomp $_;
		if ($_ =~ /^>(.*)/){
			my ($i) = split(/\s+/,$1);
			$IDs{$i}++;
		}
	}
	close FILE or die "Error closing $file: $!";
	return \%IDs;
}

