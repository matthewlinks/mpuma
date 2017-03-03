#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use lib "$ENV{'mPUMA'}/lib";	# this allows us to find the module using the lone environment variable we require to be configured by the user
use mPUMA;

# Pipeline parameters: Initializing Assembly/Analysis Processes
my ($blast_db,$working_dir,$library_file,$name_of_assembly,$VERBOSE,$PROTEIN,$down_sample,$initial_OTU_file,
$GENESPRING,$MOTHUR,$SPSN,$UNIFRAC,$R,$MEGAN,$TAXONOMICTERM,$CHIMERACHECK,$ASSEMBLER,$MAPPING

);
my @options = (
	'a=s',		\$ASSEMBLER,
	'mapping=s',	\$MAPPING,
	'c=s',		\$CHIMERACHECK,	# this can be set to 0 to skip using the chimera report information. The chimera analysis will be performed regardless of this var
	'd=s',		\$blast_db,
	'w=s',		\$working_dir,
	'n=s',		\$name_of_assembly,
	'l=s',		\$library_file,	# This would be the 2-column input which maps SFFs to libraries 
	'v=s',		\$VERBOSE,
	'downSample=s',	\$down_sample,
	'protein=s',	\$PROTEIN,
	'i=s',		\$initial_OTU_file,
	'SpSn=s',	\$SPSN,		# option to turn ON generation of Sp and Sn calculations
	'mothur=s',	\$MOTHUR,	# option to turn off generating and calling MOTHUR for stats
	'megan=s',	\$MEGAN,	# option to turn off generation of input files for MEGAN
	'genespring=s',	\$GENESPRING,	# option to turn off generation of files for GeneSpring
	'unifrac=s',	\$UNIFRAC,	# option to turn ON generation of files for UniFrac
	'r=s',		\$R,		# option to turn ON generation of files for R (library comparisons etc) 
	't=s',		\$TAXONOMICTERM	# option to set the taxonomic term level for the profiles which are generated off of the classifier / MEGAN input
);
&GetOptions(@options);
$CHIMERACHECK		= 1 unless defined $CHIMERACHECK; # default is to use the chimera checking. 
$blast_db 		= $mPUMA::cpnDB_nr_nuc unless defined $blast_db; # take the default unless the user knows better
$mPUMA::VERBOSE		= 1 if(defined($VERBOSE) eq 1);
$PROTEIN		= 1 unless defined $PROTEIN;
$down_sample		= 1 unless defined $down_sample;	# this will down sample to the minimum # reads in a library

$MOTHUR			= 1 unless defined $MOTHUR;
$MEGAN			= 1 unless (defined $MEGAN)or(defined $TAXONOMICTERM);
$GENESPRING		= 1 unless defined $GENESPRING;
$UNIFRAC		= 1 unless defined $UNIFRAC;
$TAXONOMICTERM		= 'phylum' unless defined $TAXONOMICTERM; # this needs to be set to a default value AFTER the MEGAN param is set to its default. 
$SPSN			= 0 unless defined $SPSN;
$R			= 0 unless defined $R;
$ASSEMBLER		= 'gsAssembler' unless defined $ASSEMBLER;
if (!(defined $MAPPING)){
	if ($ASSEMBLER eq 'gsAssembler'){
		$MAPPING = 'bowtie2';
	}elsif($ASSEMBLER eq 'Trinity'){
		$MAPPING = 'bowtie2';
	}else{
		die "Error in unsupported assembler mode";
	}
}

usage() unless (
	defined $library_file
	and defined $working_dir
	and defined $name_of_assembly
	and defined $ASSEMBLER
	and defined $MAPPING
);

# what is the primary output from the assembly exercise
if ($ASSEMBLER eq 'gsAssembler'){
	$initial_OTU_file	= '454Isotigs.fna' unless defined $initial_OTU_file; # this is a param to allow switching between 454Isotigs.fna and 454AllContigs.fna
}elsif ($ASSEMBLER eq 'Trinity'){
	$initial_OTU_file = 'Trinity.fasta' unless defined $initial_OTU_file;
}else{
	die "In an unsupported assembler case";
}


mPUMA::checkWorkingDir($working_dir);
die "Error $name_of_assembly does not seem to be a directory within $working_dir" unless (-d join('/',$working_dir,$name_of_assembly));

sub usage {

	die <<"USAGE";

Usage: $0 -l library_file -n AssemblyName -w working_dir

	OPTIONS:
			-a string		Default: $ASSEMBLER (this is the method for OTU consensus formation)
			-mapping string		Default: $MAPPING (this is the method for read mapping against the assembly)
			-c boolean		Default: 1 - this uses the C3 chimera checking report. From some ecological niches this may be overly stringent.
			-h /path/to/bt2base	Default: undefined - this can be set to point to a bowtie 2 base which is used to EXCLUDE OTUs which are actually off-target amplifications of the host genome
			-v boolean		Default: $mPUMA::VERBOSE
			-d /path/reference_db	Default: $mPUMA::cpnDB_nr_nuc
			-downSample integer	Default: $down_sample 1 - means downsample to the smallest library size, if you want to NOT downSample set this to a number largest than the size of your largest library (i.e. 1000000)
			-protein boolean	Default: $PROTEIN
			-i initialOTUmsf	Default: $initial_OTU_file (suggest that if you are changing it you probably want 454AllContigs.fna instead)
			-SpSn boolean		Default: $SPSN
			-mothur	boolean		Default: $MOTHUR
			-megan boolean		Default: $MEGAN
			-genespring boolean	Default: $GENESPRING
			-unifrac boolean	Default: $UNIFRAC
			-r boolean		Default: $R
			-t taxonomicTerm	Default: $TAXONOMICTERM
USAGE
}

die "\$ENV{'mPUMA'} environment variable is undefined.\nPlease add the environemnt variable mPUMA" unless defined $ENV{'mPUMA'};

# Parse the library file to figure out which SFFs correspond to which libraries
my $library_info = mPUMA::parse_library_file( $library_file);	# parse_library_info must ensure no non-word chars in the library name
my @NGS_file_list;
foreach my $lib (keys %{$library_info}){
	push (@NGS_file_list,@{$library_info->{$lib}});
}
my $topLevelDir = join('/',$working_dir,$name_of_assembly);
my $assemblyDir = join('/',$topLevelDir,'assembly');

# generate a non-redundant set of OTU
# the default value for initial_OTU_file are set by assembly method above so this should be the same for both methods...
my ($nr_OTU_fasta,$cluster_report,$noChimerasFile,$chimericIDs) = mPUMA::post_assembly_cleanup($assemblyDir,$CHIMERACHECK,$initial_OTU_file);

# BLASTX the clean non-chimeric OTUs vs cpndb nr 
my $blastx_seqclean_file = mPUMA::generateBLASTX( $noChimerasFile, 100, 100);

# Using the BLASTX results solve for the orientation and translation of each OTUs
my ($oriented_nrOTUs,$aa_oriented_nrOTUs) = mPUMA::translationOfStranded( $noChimerasFile, $blastx_seqclean_file);

# Cluster the OTU, retreive the longest OTU and make a protein FASTA of those non-redundant OTU
#my ($seqclean_aa_CDHitClustFnaFile, $seqclean_aa_CDHitClustReportFile) = mPUMA::cluster_cd_hit( $seqclean_aa_file);
my ($aa_oriented_CDHitClustFnaFile, $aa_oriented_CDHitClustReportFile) = mPUMA::cluster_cd_hit( $aa_oriented_nrOTUs);

# Running wateredBlast on 454Isotigs.fna query file and cpndb_nr database target file.
# Checks the waterBlast output file for errors.
die "Error $0 uses the waterBLAST program, which invokes BLAST. Therefore the -d \$blast_db parameter database file must exist and has to be populated with data. Cannot successfully invoke waterBlast without a wbTarget file. (database file)" 
	unless defined $blast_db;

my ($nr_OTU_fasta_WB,$status) = mPUMA::waterBlastIT($nr_OTU_fasta,$blast_db,'blastn');
if($status == 1){
	warn "\nThe waterBlast output file " . $nr_OTU_fasta_WB . " has been created successfully no errors reported....\n\nOutput File: " . $nr_OTU_fasta_WB . "\n\nProcessing with waterBlast Complete....\n\n" if ($VERBOSE);
}else{	
	die "wateredBlast had an error associated with file: " . $nr_OTU_fasta_WB . "\n";
}


# this is where read mapping should diverge based on the 2 methods

# make an ACE toc file
my $mapToc;
if ($MAPPING eq 'gsAssembler'){
	$mapToc = mPUMA::make_ace_toc($assemblyDir);
	# compress ace files if the TOC is in place and looks ok
	mPUMA::check_for_ace_compression($assemblyDir);
}elsif($MAPPING eq 'bowtie2'){
	# if gsAssembler was used but we want to use bowtie 2 for mapping then we have some prep to do
	# if this was assembled by Trinity then this should recognize there is no work to do and just return
 	prepFastqBasedAssembly($topLevelDir,@NGS_file_list);
	$mapToc = mPUMA::make_sam_toc($topLevelDir,$assemblyDir,$nr_OTU_fasta);
}else{
	die "Error in an unsupported form of assembler";
}

# Creates the total read set file "total.fna" from sffinfo (reads information from the SFF file).
# # Checks the format of the total readset file "total.fna" output file for errors.
my $totalReadSet;
if ($ASSEMBLER eq 'gsAssembler'){
	($totalReadSet, $status) = mPUMA::sffinfoFna( $assemblyDir,@NGS_file_list);
	if($status eq 1){
		warn  "\nThe total read set file \"" . $totalReadSet . "\" has been created successfully no errors reported....\n\n" if ($VERBOSE);
	}
	else{
		die "sffinfoFna had an error associated with file: " . $totalReadSet . ": $!\n";
	}
}elsif($ASSEMBLER eq 'Trinity'){
	$totalReadSet = join('/',$assemblyDir,'total.fna');
	die "Error lost $totalReadSet" unless (-e $totalReadSet);
}else{
	die "Error you are in an unsupported assembler case";
}

# create a toc file which is the breakdown of where every read ended up in the nucleotide assembly
# this should include the SeqCleaning and the CD-Hitting 
# READ,FINAL_OTU_ID,RATIONALE
# 
my $totalToc = mPUMA::make_total_toc($assemblyDir,$totalReadSet,$mapToc,$noChimerasFile,$cluster_report,$chimericIDs,$CHIMERACHECK);

die "Error total TOC file appears to be messed up: $totalToc" unless (-s $totalToc);


# Running wateredBlast on total.fna query file and cpndb_nr database target file. 
# Checks the waterBlast output file for errors.
my $totalReadSetWB;
($totalReadSetWB, $status) = mPUMA::waterBlastIT( $totalReadSet, $blast_db, 'blastn');
if($status == 1){
	warn "\nThe waterBlast output file [" . $totalReadSetWB . "] has been created successfully no errors reported....\n\nOutput File: " . $totalReadSetWB . "\n\nProcessing with waterBlast Complete....\n\n" if ($VERBOSE);
}else{
	die "wateredBlast had an error associated with file: " . $totalReadSetWB . ": $!\n";
}

# Call the Classifier on both the total set of READs and the assembly
my $totalReadSetClassifier = mPUMA::run_classifier_on_fasta($totalReadSet);
my $nr_OTU_fasta_classifier = mPUMA::run_classifier_on_fasta($nr_OTU_fasta);

# Creates the nucleotide / protien and project directory if they don't exist.
my $nuc_dir = join('/', $working_dir, 'nucleotide');
unless(-d $nuc_dir){
	mkdir($nuc_dir, 0777) or die "Can't make directory $nuc_dir: $!";
}
my $library_dir_nuc = mPUMA::createLibraries( $nuc_dir, $library_info);

my ($library_dir_pep,$pep_dir);
if ($PROTEIN){
	$pep_dir = join('/', $working_dir, 'protein');
	unless(-d $pep_dir){
		mkdir($pep_dir, 0777) or die "Can't make directory $pep_dir: $!";
	}
	$library_dir_pep = mPUMA::createLibraries( $pep_dir, $library_info);
}

# Creates "reads.ids" file in each project directory which lists the read ids specific to that library.
# Checks the "reads.ids" file in each project directory for errors and computes the median number of 
# calculates the median number of reads / library
warn "\nGenerate Read IDs files for conversion from assembled data into diversity data for analysis....\n\n" if ($VERBOSE);
my @readCountsNuc; # array to hold the numbr of reads / library
if ((require Parallel::Loops)and($mPUMA::cpu_cores)){
        # do the wateredBLAST of total.fna in parallel
        my $parallel = Parallel::Loops->new($mPUMA::cpu_cores);
        my %readCounts;
        $parallel->share(\%readCounts); # make sure that these are visible in the children --- this might be creating this as a semaphore but I'm not sure without looking this up in the perl code

        # run these in parallel
	my @jobs = keys %{$library_info};
        $parallel->foreach(\@jobs,sub {
		my $lib = $_;
		my $lDir = join('/',$library_dir_nuc,$lib);
#                my $numReads = mPUMA::generateReadIds($lDir,$library_info->{$lib});
                my $numReads;
		if ($MAPPING eq 'gsAssembler'){
                	$numReads = mPUMA::generateReadIds($lDir,$library_info->{$lib});
		}elsif($MAPPING eq 'bowtie2'){
                	$numReads = mPUMA::generateReadIdsFastq($topLevelDir,$lDir,$library_info->{$lib});
		}else{
			die "In unsupported read mapping";
		}
		$readCounts{$lib} = $numReads;
        });
        foreach my $lib (keys %readCounts){
                push(@readCountsNuc,$readCounts{$lib});
        }
        undef $parallel;
}else{
	foreach my $lib (keys %{$library_info}){
		my $lDir = join('/',$library_dir_nuc,$lib);
#		my $numReads = mPUMA::generateReadIds( $lDir,$library_info->{$lib});
                my $numReads;
		if ($MAPPING eq 'gsAssembler'){
                	$numReads = mPUMA::generateReadIds($lDir,$library_info->{$lib});
		}elsif($MAPPING eq 'bowtie2'){
                	$numReads = mPUMA::generateReadIdsFastq($topLevelDir,$lDir,$library_info->{$lib});
		}else{
			die "In unsupported read mapping";
		}
		push(@readCountsNuc,$numReads);
	}
}

my $medianNumReadsNuc = mPUMA::find_median( @readCountsNuc);
warn "The median number of Reads for nucleotide is $medianNumReadsNuc\n\n" if ($VERBOSE);

my $minNumReadsNuc = mPUMA::find_minimum(@readCountsNuc);
warn "The min number of Reads in a library is $minNumReadsNuc" if $VERBOSE;

# down sample the nucleotide data
if ($down_sample > 0){
	foreach my $lib (keys %{$library_info}){
		my $lDir = join('/',$library_dir_nuc,$lib);
		if ($down_sample == 1){
			mPUMA::subsample_library($lDir,$minNumReadsNuc);
		}else{
			mPUMA::subsample_library($lDir,$down_sample);
		}
	}
}

my @readCountsPep; # array to hold the numbr of reads / library
if($PROTEIN){
	if ((require Parallel::Loops)and($mPUMA::cpu_cores)){
	        # do the wateredBLAST of total.fna in parallel
	        my $parallel = Parallel::Loops->new($mPUMA::cpu_cores);
	        my %readCounts;
	        $parallel->share(\%readCounts); # make sure that these are visible in the children --- this might be creating this as a semaphore but I'm not sure without looking this up in the perl code
	
	        # run these in parallel
	        my @jobs = keys %{$library_info};
	        $parallel->foreach(\@jobs,sub {
			my $lib = $_;
	 		my $pepDir = join('/',$library_dir_pep,$lib);
	 		my $nucDir = join('/',$library_dir_nuc,$lib);

			my $nucReads = join('/',$nucDir,$mPUMA::read_ids_file);
			my $pepDest =  join('/',$pepDir,$mPUMA::read_ids_file);
	
			# only symlink if the protein data isn't there
			unless (-e $pepDest){
				symlink($nucReads,$pepDest) or die "Error creating symlink $nucReads,$pepDest: $!";
			}

			# if we are down sampling then ensure that those links exist too
			if ($down_sample > 0){
				$nucReads = join('/',$nucDir,$mPUMA::read_sub_file);
				$pepDest =  join('/',$pepDir,$mPUMA::read_sub_file);
		
				# only symlink if the protein data isn't there
				unless (-e $pepDest){
					symlink($nucReads,$pepDest) or die "Error creating symlink $nucReads,$pepDest: $!";
				}
			}
	        });
	        foreach my $lib (keys %readCounts){
	                push(@readCountsPep,$readCounts{$lib});
	        }
	        undef $parallel;
	}else{
		foreach my $lib (keys %{$library_info}){
			my $pepDir = join('/',$library_dir_pep,$lib);
                        my $nucDir = join('/',$library_dir_nuc,$lib);

                        my $nucReads = join('/',$nucDir,$mPUMA::read_ids_file);
                        my $pepDest =  join('/',$pepDir,$mPUMA::read_ids_file);
        
                        # only symlink if the protein data isn't there
                        unless (-e $pepDest){
                                symlink($nucReads,$pepDest) or die "Error creating symlink $nucReads,$pepDest: $!";
                        }

                        # if we are down sampling then ensure that those links exist too
                        if ($down_sample > 0){
                                $nucReads = join('/',$nucDir,$mPUMA::read_sub_file);
                                $pepDest =  join('/',$pepDir,$mPUMA::read_sub_file);

                                # only symlink if the protein data isn't there
                                unless (-e $pepDest){
                                        symlink($nucReads,$pepDest) or die "Error creating symlink $nucReads,$pepDest: $!";
                                }
                        }


			my $lDir = join('/',$library_dir_pep,$lib);
			my $numReads = mPUMA::generateReadIds( $lDir,$library_info->{$lib});
			push(@readCountsPep,$numReads);
		}
	}
}


# Create diversity files for each MID directory. Scaled to the median number of reads in the library.
my @diversityFilesNuc;	# the list of files
my %lib2DiversityFilesNuc;	# a hash to get the file for a specific library
warn "Generating the nucleotide-based diversity files....\n" if($VERBOSE);
if ((require Parallel::Loops)and($mPUMA::cpu_cores)){
        # do the wateredBLAST of total.fna in parallel
        my $parallel = Parallel::Loops->new($mPUMA::cpu_cores);
	my %files;
	$parallel->share(\%files);
	my @jobs = keys %{$library_info};
	$parallel->foreach(\@jobs,sub {
		my $lib = $_;
		my $lDir = join('/',$library_dir_nuc,$lib);
		$files{$lib} = mPUMA::generateDiversityToc($lDir,$totalToc,$nr_OTU_fasta_WB);
	});
	foreach my $l (keys %files){
		$lib2DiversityFilesNuc{$l} = $files{$l};
                push(@diversityFilesNuc,$files{$l});
	}
	undef $parallel;
}else{
	foreach my $lib (keys %{$library_info}){
		my $lDir = join('/',$library_dir_nuc,$lib);
		my $divFile = mPUMA::generateDiversityToc($lDir,$totalToc,$nr_OTU_fasta_WB);
		$lib2DiversityFilesNuc{$lib} = $divFile;
		push(@diversityFilesNuc,$divFile);
	}
}

# Create diversity files for each MID directory. Scaled to the median number of reads in the library.
my @diversityFilesPep;	# the list of files
my %lib2DiversityFilesPep;	# a hash to get the file for a specific library
if ($PROTEIN){
	warn "Generating the protein-based diversity files....\n" if($VERBOSE);
	if ((require Parallel::Loops)and($mPUMA::cpu_cores)){
	        # do the wateredBLAST of total.fna in parallel
	        my $parallel = Parallel::Loops->new($mPUMA::cpu_cores);
	        my %files;
	        $parallel->share(\%files);
	        my @jobs = keys %{$library_info};
	        $parallel->foreach(\@jobs,sub {
			my $lib = $_;
	                my $lDir = join('/',$library_dir_pep,$lib);
			$files{$lib} = mPUMA::generateDiversityToc($lDir,$totalToc,$nr_OTU_fasta_WB,$aa_oriented_CDHitClustReportFile);
	        });
	        foreach my $l (keys %files){
			push(@diversityFilesPep,$files{$l});
			$lib2DiversityFilesPep{$l} = $files{$l};
	        }
	        undef $parallel;
	}else{
		foreach my $lib (keys %{$library_info}){
			my $lDir = join('/',$library_dir_pep,$lib);
			my $divFile = mPUMA::generateDiversityToc($lDir,$totalToc,$nr_OTU_fasta_WB,$aa_oriented_CDHitClustReportFile);
			$lib2DiversityFilesPep{$lib} = $divFile;
			push(@diversityFilesPep,$divFile);
		}
	}
}

### ADD step for mothur input file
my @libs;
foreach my $lib (keys %{$library_info}){
	push(@libs,$lib);
}
if ($MOTHUR){
	warn "\nGenerating the nucleotide MOTHUR input file....\n\n" if($VERBOSE);
	my ($mothur_nuc_file, $mothur_nuc_files_dir) = mPUMA::generateInputMOTHUR($nuc_dir,$library_dir_nuc,\@libs);
	
	warn "\nGenerating nucleotide output files for MOTHUR....\n\n" if($VERBOSE);
	mPUMA::processInputMOTHUR( $mothur_nuc_file, $mothur_nuc_files_dir);

	if ($PROTEIN){
		warn "\nGenerating the protein MOTHUR input file....\n\n" if($VERBOSE);
		my ($mothur_pep_file, $mothur_pep_files_dir);
		($mothur_pep_file, $mothur_pep_files_dir) = mPUMA::generateInputMOTHUR( $pep_dir, $library_dir_pep,\@libs);
		
		warn "\nGenerating protein output files for MOTHUR....\n\n" if($VERBOSE);
		mPUMA::processInputMOTHUR( $mothur_pep_file, $mothur_pep_files_dir);
	}
}

# Make Megan output which is based off of the Classifier results and the abundance in the diversity files
if ($MEGAN){
	warn "\nGenerating the nucleotide MEGAN Classifier files....\n\n" if($VERBOSE);
	my %meganNucFiles;
	if ((require Parallel::Loops)and($mPUMA::cpu_cores)){
	        # do the wateredBLAST of total.fna in parallel
	        my $parallel = Parallel::Loops->new($mPUMA::cpu_cores);
	        $parallel->share(\%meganNucFiles);
	        my @jobs = keys %{$library_info};
	        $parallel->foreach(\@jobs,sub {
			my $lib = $_;
			$meganNucFiles{$lib} = mPUMA::generateClassifierInputMEGAN( $nuc_dir, $lib, $lib2DiversityFilesNuc{$lib}, $nr_OTU_fasta_classifier);
			mPUMA::generateClassifierProfilesFromMegan($nuc_dir, $lib,$meganNucFiles{$lib},$TAXONOMICTERM) if (defined $TAXONOMICTERM);
		});
	        undef $parallel;
	}else{
		foreach my $lib (keys %{$library_info}){
			$meganNucFiles{$lib} = mPUMA::generateClassifierInputMEGAN( $nuc_dir, $lib, $lib2DiversityFilesNuc{$lib}, $nr_OTU_fasta_classifier);
			mPUMA::generateClassifierProfilesFromMegan($nuc_dir, $lib,$meganNucFiles{$lib},$TAXONOMICTERM) if (defined $TAXONOMICTERM);
		}
	}

	if($PROTEIN){
		# Generate the Megan CSV files based on the protein space
		warn "\nGenerating the protein MEGAN Classifier files....\n\n" if($VERBOSE);
		my %meganPepFiles;
		if ((require Parallel::Loops)and($mPUMA::cpu_cores)){
		        # do the wateredBLAST of total.fna in parallel
		        my $parallel = Parallel::Loops->new($mPUMA::cpu_cores);
		        $parallel->share(\%meganPepFiles);
		        my @jobs = keys %{$library_info};
		        $parallel->foreach(\@jobs,sub {
				my $lib = $_;
				$meganPepFiles{$lib} = mPUMA::generateClassifierInputMEGAN( $pep_dir, $lib, $lib2DiversityFilesPep{$lib}, $nr_OTU_fasta_classifier);
				mPUMA::generateClassifierProfilesFromMegan($pep_dir, $lib,$meganPepFiles{$lib},$TAXONOMICTERM) if (defined $TAXONOMICTERM);
			});
		        undef $parallel;
		}else{
			foreach my $lib (keys %{$library_info}){
				$meganPepFiles{$lib} = mPUMA::generateClassifierInputMEGAN( $pep_dir, $lib, $lib2DiversityFilesPep{$lib}, $nr_OTU_fasta_classifier);
				mPUMA::generateClassifierProfilesFromMegan($pep_dir, $lib,$meganPepFiles{$lib},$TAXONOMICTERM) if (defined $TAXONOMICTERM);
			}
		}
	}
}

if ($GENESPRING){
	
	# To load data into GeneSpring you need to have diversity files which have a row for each and every OTU even if it wasn't seen in that library

	# Create a list of isotigs that exist in the diversity file
	my $OTU_nuc_list_file = mPUMA::generateIsotigList( $nuc_dir,\@diversityFilesNuc);
	
	# Creates template technology files for OTU extraction used by the convertTech sub routine.
	my $technology_nuc_file = mPUMA::generateTechFile( $OTU_nuc_list_file, 10000);
	
	# Converts the template files into technology (*.tsv) files usable by gene_spring.
	# Checks the convertedTech output formated file for errors.
	my $genespring_nuc_dir;
	($genespring_nuc_dir,$status) = mPUMA::convertTech( $nuc_dir, $technology_nuc_file,$library_dir_nuc,\@libs);
	# $nuc_dir,$library_dir_nuc,\@libs
	if($status eq 1){
		warn "\nAll \"convertedTSV\" files were successfully created no errors reported....\n\n" if ($VERBOSE);
	} 
	else{
		die "convertTech had an error $nuc_dir, $technology_nuc_file (",join(',',@diversityFilesNuc),")\n$status";
	}
	

	if ($PROTEIN){
		my $genespring_pep_dir;
		my $OTU_pep_list_file = mPUMA::generateIsotigList( $pep_dir,\@diversityFilesPep) if ($PROTEIN);
		my $technology_pep_file = mPUMA::generateTechFile( $OTU_pep_list_file, 10000) if ($PROTEIN);
	
		($genespring_pep_dir,$status) = mPUMA::convertTech( $pep_dir, $technology_pep_file,$library_dir_pep,\@libs);
		
		if($status eq 1){
			warn "\nAll \"convertedTSV\" files were successfully created no errors reported....\n\n" if ($VERBOSE);
		} 
		else{
			die "convertTech had an error\n$status";
		}
	}	
}

# generate input files for UniFrac
if ($UNIFRAC){
	warn "\nGenerating the nucleotide UNIFRAC input file....\n\n" if($VERBOSE);
	my $nucUnifracDir = mPUMA::generateInputUNIFRAC( $nuc_dir, \@diversityFilesNuc,$oriented_nrOTUs,'Nucleotide');
	if ($PROTEIN){
		warn "\nGenerating protein the UNIFRAC input file....\n\n" if($VERBOSE);
		my $pepUnifracDir = mPUMA::generateInputUNIFRAC( $pep_dir, \@diversityFilesPep,$aa_oriented_CDHitClustFnaFile,'Protein');
	}
}

# Computes the Sensitivity and Specificity metrics using the total.fna file and cpndb_nr database file. 
if ($SPSN){
	# Checks the "Specificity_Sensitivity" output file for errors.
	
	# This is currently commented out for performance 
	# As this is completely independent it should be forked off
	
	my ($snspInputFile, $snspOutputFile);
	($snspInputFile, $snspOutputFile, $status) = mPUMA::calcSpSnMetrics( $assemblyDir, $totalReadSet, $nr_OTU_fasta,$totalReadSetWB,$nr_OTU_fasta_WB);
	
	if($status eq 1){
		warn "\nThe validation metrics file " . $snspInputFile . " has been created successfully no errors reported....\n\nOutput File: " . $snspInputFile . "\n\n" if ($VERBOSE);
	}
	else{
		die "calcSpSnMetrics had an error associated with file: " . $snspInputFile . ": $!\n";
	}
}



### done
exit 0;
