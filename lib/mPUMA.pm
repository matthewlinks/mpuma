#############################################################
# mPUMA.pm - Microbial Profiling Using Metagenomic Assembly #
# Metagenomic assembly and analysis of CPN60 454 data       #
#############################################################

package mPUMA;
require Exporter;
use strict;
use IPC::Open2;
use POSIX;
use File::Basename;
use File::Find;
use IO::Uncompress::Gunzip;
use Bio::SeqIO;
use Cwd qw(chdir abs_path);
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);
# Exports all Global Parameters and Commands.
@EXPORT = qw($newAssembly $addRun $runProject $formatdb $seqclean $cdbfasta $cdbyank $cdHitEST $blastall $mothur $watered_blast $calculate_sp_sn_metrics $parseCDHitEst $convertTechCmd $sffinfo $fixStrand $generateDiversity $createMeganInput $createMeganBlastxInput $createMothurInput $cpu_cores $cpu_seqclean $min_length $min_identity $primer_file $cpnDB_nr_nuc $cpnDB_nr_pep $num_descs $read_ids_file $VERBOSE parse_library_file newAssemb addRun runProject runProjectCheck runProjectError runAssembly runTrinity post_assembly_cleanup seqcleanFile cluster_cd_hit waterBlastIT waterBlastCheck sffinfoFna sffsToFna sffinfoFnaCheck calcSpSnMetrics calcSpSnMetricsCheck  solveOrientation createLibraries generateReadsIds find_median generateDiversityCDHit getNUM_CPU generateDiversityOnlyIsotigs generateIsotigList generateTechFile convertTech convertTechCheck generateBLASTX generateFastaMEGAN generateBlastxMEGAN generateCSVInputMEGAN processCSVInputMEGAN generateInputMOTHUR generateInputUNIFRAC processInputMOTHUR split_fasta_n_chunks combineWateredBlastFiles nameWaterBlastResultFile checkWorkingDir prepFastqBasedAssembly );


# Globals
our ($newAssembly, $addRun, $runProject, $trinity, $formatdb, $seqclean, $cdbfasta, $cdbyank, $cdHitEST, $cdHit, $blastall, $mothur, $translateStranded, $watered_blast, $calculate_sp_sn_metrics, $parseCDHitEst, $convertTechCmd, $sffinfo, $fixStrand, $generateDiversity, $createMeganInput, $createMeganBlastxInput, $createMeganCSVInput, $createMothurInput, $createUnifracInput,$makeToc, $createMeganClassifierInput,$fastTree, $tCoffee, $classifier, $bowtie2, $bowtie2_build

);

# Edit these paths to ensure that mPUMA can find the executables it needs
$newAssembly			= '/usr/local/rig/bin/newAssembly';
$addRun				= '/usr/local/rig/bin/addRun';
$sffinfo			= '/usr/local/rig/bin/sffinfo';
$runProject			= '/usr/local/rig/bin/runProject';
$trinity			= '/usr/local/trinity/Trinity.pl';
$formatdb			= '/usr/local/blast/bin/formatdb';
$seqclean			= '/usr/local/seqclean/seqclean';
$cdbfasta			= 'cdbfasta';
$cdbyank			= 'cdbyank';
$cdHitEST			= '/usr/local/cd-hit/cd-hit-est';
$cdHit				= '/usr/local/cd-hit/cd-hit';
$blastall			= '/usr/local/blast/bin/blastall';
$mothur				= '/usr/local/bin/mothur';
$tCoffee			= 't_coffee';
$fastTree			= "FastTreeMP";			# http://meta.microbesonline.org/fasttree/ install the multi-threaded version
$classifier			= "$ENV{mPUMA}/classifier/rdp_classifier_2.2/rdp_classifier-2.2.jar";
$bowtie2			= '/usr/local/bowtie2/bowtie2';
$bowtie2_build			= '/usr/local/bowtie2/bowtie2-build';

# Global Parameters
our ($cpu_cores, $cpu_seqclean, $min_length, $min_identity, $primer_file, $cpnDB_nr_nuc, $cpnDB_nr_pep, $num_descs, $read_ids_file, $VERBOSE,
$read_sub_file,$diversity_file
);

# edit these parameters if you want to affect the # cpus mPUMA will use or the assembly parameters
$cpu_cores	= 30;	# -c parameter to newbler
$min_length	= 137;	# -ml parameter to newbler
$min_identity	= 90;	# -mi parameter to newbler
$cpu_seqclean	= 16; 	# must be <= 16 due to limit within seqclean
#Classic
$primer_file    = "$ENV{mPUMA}/reference_fasta/primers-enumerated-3primeN.msf";
# rpoB
#$primer_file	= "$ENV{mPUMA}/reference_fasta/rpoB-primers.fna.enumerated";
$cpnDB_nr_nuc   = "$ENV{mPUMA}/reference_fasta/cpndb_nr_20140327";
$cpnDB_nr_pep   = "$ENV{mPUMA}/reference_fasta/cpndb_nr_pep_20140327";
$VERBOSE	= 0;

# These are internal programs within mPUMA and you should not need to edit these unless you are a developer
$translateStranded		= "$ENV{mPUMA}/bin/translate_stranded.pl";
$watered_blast			= "$ENV{mPUMA}/bin/watered_blast.pl";
$calculate_sp_sn_metrics	= "$ENV{mPUMA}/bin/calculate_Sp_Sn_metrics.pl";
$parseCDHitEst			= "$ENV{mPUMA}/bin/parse_cd-hit-est.pl";
$convertTechCmd 		= "$ENV{mPUMA}/bin/convert_technology.pl";
$fixStrand			= "$ENV{mPUMA}/bin/fix_strandedness.pl";
$generateDiversity		= "$ENV{mPUMA}/bin/generate_diversity_from_toc.pl";
$createMeganInput		= "$ENV{mPUMA}/bin/create_input_megan.pl";
$createMeganBlastxInput 	= "$ENV{mPUMA}/bin/create_blastx_input_megan.pl";
$createMeganCSVInput		= "$ENV{mPUMA}/bin/create_CSV_input_megan.pl";
$createMeganClassifierInput	= "$ENV{mPUMA}/bin/create_input_megan_classifier.pl";
$createMothurInput		= "$ENV{mPUMA}/bin/create_input_mothur.pl";
$createUnifracInput		= "$ENV{mPUMA}/bin/create_input_unifrac.pl";
$makeToc			= "$ENV{mPUMA}/bin/make_toc.pl";
$num_descs			= 5;	# used in wateredBLAST searches
$read_ids_file			= 'reads.ids';	
$read_sub_file			= 'subsample.ids';
$diversity_file			= 'diversity';

sub checkWorkingDir {
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	die "Error $dir does not look like a fully qualified path. Try specifying it from the root of the filesystem" unless($dir =~ /^[\/\~]/);
}

# sub routine to split a fasta file into pieces
# number of pieces can be passed as an argument or defaults to $mPUMA::cpu_cores
sub split_fasta_n_chunks {
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $num_pieces = shift;
	$num_pieces = $cpu_cores unless defined $num_pieces;

	my ($filename,$base_dirs) = fileparse($file);
	my $prefix = $filename . '.split';

	my @files;
	my $number_non_zero_files = 0;
	for(my $i = 1; $i <= $num_pieces; $i++){
		my $f = $file . '.split.' . $i;
		push(@files,$f);
		$number_non_zero_files++ if (-s $f);
	}

	if ($number_non_zero_files == $num_pieces){
		warn "It appears that the file split was already done and is assumed to be correct / accurrate";
	}else{
		system("$ENV{mPUMA}/bin/split_fasta.pl",
			'-i',	$file,
			'-c',	$num_pieces,
			'-p',	$prefix,
			'-o',	$base_dirs) == 0 or die "Error splitting $file into $num_pieces pieces: $!";
	}
	return \@files;
}

# sub routine to combine multiple watered BLAST reports together
sub combineWateredBlastFiles {
	my $files = shift;
	die "Error lost files to combine" unless defined $files;
	my $wbQuery = shift;
	die "Error need the query file" unless defined $wbQuery;
	my $wbDB = shift;
	die "Error need the query file" unless defined $wbDB;

	my $wbOutputFile = nameWaterBlastResultFile($wbQuery,$wbDB);
	if (-s $wbOutputFile){
		warn "Error $wbOutputFile exists already so not going to replace it";
		return $wbOutputFile;
	}
	open(OUTPUT,">$wbOutputFile") or die "Error opening $wbOutputFile for writting: $!";
	my $written_header = 0;
	foreach my $file (@$files){
		open(INPUT,"<$file") or die "Error opening $file for reading: $!";
		my $header = <INPUT>; 
		chomp $header;
		print OUTPUT $header,"\n" unless $written_header;
		$written_header++;
		while(<INPUT>){
			print OUTPUT $_;
		}
		close INPUT or die "Error closing $file: $!";
	}
	close OUTPUT or die "Error closing $wbOutputFile: $!";
	return $wbOutputFile;
}

#############################################################################
# parse_library_file - Simple parser which reads in a 2-column file which   #
# defines which SFFs are for which libraries and is tab delimited           #
#                                                                           #
# E.g.                                                                      #
#                                                                           #
# LibA	/path/to/SFF1.sff                                                   #
# LibB	/path/to/SFF2.sff                                                   #
# LibB	/path/to/SFF3.sff                                                   #
#############################################################################
# Input Parameters:                                                         #
#                   $libfile  - The library file for parsing SFFs           #
# Output Parameters:                                                        #
#                   \%info - Hash reference pointing to the SFF filename    #
#                            paths                                          #
#                                                                           #
#############################################################################

sub parse_library_file {
	my $libfile = shift;
	die "Error lost library file" unless defined $libfile;
	open(FILE,"<$libfile") or die "Error opening $libfile: $!";
	my %info;
	while(<FILE>){
		chomp $_;
		$_ =~ s/\r//g; # try and strip characters which could arise from old text editors
		my ($lib,$ngsFile) = split(/\t/,$_);
		if ($lib =~ /^[\w\_\-]+$/){
			push(@{$info{$lib}},$ngsFile);
		}else{
#			die "Error lib ($lib)does not match only [\\w\\_\\-]+: $_";
			warn "Error lib ($lib)does not match only [\\w\\_\\-]+: $_";
			push(@{$info{$lib}},$ngsFile);
		}
	}
	close FILE or die "Error closing $libfile: $!";
	return \%info;
}

#############################################################################
# newAssembly - Creates the assembly directory and populates it with files  #
# from the cdna library.                                                    #
#############################################################################
# Input Parameters:                                                         #
#                   $newblerDir - The Newbler project directory             #
# Output Parameters:                                                        #
#                   $assemblyDir- The assembly directory                    #
#                                                                           #
#############################################################################

sub newAssembly{
	#newAssembly -cdna /some/dir
	my $newblerDir = shift;
        die "Error lost the newbler directory" unless defined $newblerDir;

	my $assemblyDir = join('/',$newblerDir,'assembly');
	warn "\nInitializing Genome Assembly Pipeline....\n\nCreating newbler assembly directory....\n" if ($VERBOSE);

	# Creates the diversity directory if one does not exist. 
	# If it exists then exit out of sub-routine
	if(-d $assemblyDir) {
		runProjectCheck($assemblyDir);
		warn "The assembly directory exists!....\n\n" if ($VERBOSE);
		return $assemblyDir;
	}
	else{
		system($newAssembly, '-cdna', $newblerDir) == 0 or die "Error calling newAssembly –cdna $newblerDir: $?";
	}
	return $assemblyDir;
}

#############################################################################
# addRun - Adds sequence runs from the corresponding SFF files.             #
#############################################################################
# Input Parameters:                                                         #
#                   $newblerDir - The Newbler project directory             #
#                                                                           #
#############################################################################

sub addRun{
	#addRun /some/dir /some/sff/file.sff
	my $newblerDir = shift;
        die "Error lost the newbler directory" unless defined $newblerDir;
	warn "Adding SFF files to the assembly....\n" if ($VERBOSE);

	if(scalar(@_) ne 0){
		foreach my $arg (@_){
        		warn join(' ', $addRun, $newblerDir, $arg), "\n" if ($VERBOSE);
		}
        	system($addRun, $newblerDir, @_) == 0 or die "Error calling addRun $newblerDir @_: $?";
	}
	elsif(scalar(@_) eq 0){
		warn "SFF file path parameters empty!\n\n" if ($VERBOSE);
		return;
	}
}

#############################################################################
# runProject - Initiates the project assembly using newbler's runProject    #
# command.                                                                  #
#############################################################################
# Input Parameters:                                                         #
#                   $assemblyDir- The assembly directory                    #
#                                                                           #
#############################################################################

sub runProject{
	my $assemblyDir = shift;
        die "Error lost the assembly directory" unless defined $assemblyDir;

	my $Newbler454ProgFile = join('/',$assemblyDir,'454NewblerProgress.txt');
	my $isotigsLayoutFile = join('/',$assemblyDir,'454IsotigsLayout.txt');
	my $aceDIR = join('/',$assemblyDir,'ace');

	if(-s $Newbler454ProgFile and -s $isotigsLayoutFile and -d $aceDIR){
		# Checks the outcome of the project assembly for any errors and troubleshoots any issues associated 
		# with the error.
		runProjectCheck($assemblyDir);
	}
	else{
#		my $runProjectCmd = "runProject -cpu $cpu_cores -rip -acedir -m -ml $min_length -mi $min_identity -vt $primer_file";
		my $runProjectCmd = "runProject -cpu $cpu_cores -rip -acedir -m -ml $min_length -mi $min_identity";
		warn $runProjectCmd . "\n" if ($VERBOSE);
		if(-d $assemblyDir){
			chdir($assemblyDir) or die "Can't change current directory: $!";
		}
		system($runProject,
			'-cpu',	$cpu_cores,	
			'-rip',    # this doesn't actually guarantee that a read only ends up in one place in -cdna mode
			'-acedir', # Need ace files created 
			'-m',
			'-ml',	$min_length,
			'-mi',	$min_identity, 
#			'-vt',	$primer_file, # has been disabled because as best we can tell gsAssembler is not doing anything and seqclean does a great job
			'-ig',	10000,
			'-it',	10000,
			'-icc',	200
		) == 0 or die "Error calling $runProjectCmd: $?";
 
		# Checks the outcome of the project assembly for any errors and troubleshoots any issues associated 
		# with the error.
		runProjectCheck($assemblyDir);
	}

}

#############################################################################
# runProjectCheck - Checks the outcome of the project assembly for any      #
# errors and troubleshoots any issues associated with the error.            #
#############################################################################
# Input Parameters:                                                         #
#                   $assemblyDir - The assembly directory                   #
#                                                                           #
#############################################################################

sub runProjectCheck{
	my $assemblyDir = shift;
        die "Error lost the assembly directory" unless defined $assemblyDir;

	my $Newbler454ProgFile = join('/',$assemblyDir,'454NewblerProgress.txt');
	my $pipe;
	open $pipe, "-|", "/usr/bin/tail", "-n 1", $Newbler454ProgFile
	    or die "could not start tail on 454NewblerProgress.txt: $!";

	while(<$pipe>){
		if($_ =~ m/Assembly computation succeeded at: \w+ \w+ \d+ \d+:\d+:\d+ \d+/){
			warn "Assembly was Successful....\n454NewblerProgress\n" . $_ . "\n" if ($VERBOSE);
		}
		else{
			close $pipe or die $!;
			open $pipe, "-|", "/usr/bin/tail", "-n 10", $Newbler454ProgFile
	    		or die "could not start tail on 454NewblerProgress.txt: $!";
			my $errorMsg = "\nAssembly failed check $assemblyDir....\nYou may need to delete the directory, which would delete any post assembly analysis!\n\n";
			while(<$pipe>){
				$errorMsg = $errorMsg . $_;
			}
			close $pipe or die $!;
			die $errorMsg. " $!\n";# "Assembly failed check $assemblyDir....\nYou may need to delete the directory, which would delete any post assembly analysis!\n" . @error . " $!\n";
		}
	}
	close $pipe or die $!;

	my $isotigsLayoutFile = join('/',$assemblyDir,'454IsotigsLayout.txt');

	open $pipe, "-|", "/usr/bin/tail", "-n 2", $isotigsLayoutFile
	    or die "could not start tail on 454IsotigsLayout.txt: $!";

	while(<$pipe>){
 		if($_ =~ m/^#isotig/){
			warn "454IsotigsLayout\n" . $_ if($VERBOSE);
			if($_ =~ m/\#ig_thresh/ or $_ =~ m/\#it_thresh/){
				runProjectError();
			}
		}
		else{
			warn $_ if ($VERBOSE);
		}
	}
	close $pipe or die $!;

	my @listing = ();
	my $aceDIR = join('/',$assemblyDir,'ace');
	opendir(DIR, $aceDIR) or die "There is no ace directory $aceDIR: $!";
	while(my $aceFiles = readdir(DIR)){
		if ($aceFiles =~ m/.+.ace/){
			push(@listing, $aceFiles);
		}
	}
	closedir(DIR);

	if (scalar(@listing) ne 0){
		warn "\nThe ace directory and ace files have been successfully created....\n\n" if ($VERBOSE);
	}
	else{
		die "\nThe ace files have not been created successfully....\n\n";
	}
}

#############################################################################
# runProjectError - Reports an error associated with the runProject         #
#	            assembly.                                               #
#############################################################################

sub runProjectError{

die <<"ERROR";

Re-run the runProject command using the following;

runProject -cpu $cpu_cores -rip -acedir -m -ml $min_length -mi $min_identity -ig 10000 -it 10000 -icc 200 -vt $primer_file

ERROR

}

#############################################################################
# runAssembly - Initiates the project assembly using newbler's newAssembly, #
# addRun, and runProject commands.                                          #
#############################################################################
# Input Parameters:                                                         #
#                   $newblerDir - The Newbler project directory             #
# Output Parameters:                                                        #
#                   $assemblyDir - The assembly directory                   #
#                                                                           #
#############################################################################

sub runAssembly{
	
	my $newblerDir = shift;
        die "Error lost the newbler directory" unless defined $newblerDir;

	# Creates the assembly directory and populates it with files from the cdna library.
	my $assemblyDir = newAssembly($newblerDir);
	# Adds the sequence runs from the SFF files.
	addRun($newblerDir, @_);
	
	# Initiates the project assembly and records the time needed to complete the assembly process.
	# Checks the outcome of the project assembly for any errors and troubleshoots any issues associated 
	# with the error.
	runProject($assemblyDir);
	
	return $assemblyDir;
}


# simple procedure to create a Table of Contents file from all of the ACE files which can be used to identify which reads ended up in which OTUs 

sub make_ace_toc {
	my $assemblyDir = shift;
	die "Error lost assemblyDir" unless defined $assemblyDir;
	my $aceDir = join('/',$assemblyDir,'ace');
	my $aceToc = join('/',$assemblyDir,'ace.toc');
	if (!(-e $aceToc)){
		system($makeToc,
			'-d',	$aceDir,
			'-o',	$aceToc,
			'-n',	$cpu_cores) == 0 or die "Error calling $makeToc: $!";
		die "Error makeing $aceToc: $!" unless (-e $aceToc);
	}
	return $aceToc;
}

sub get_compressed_fastq_folder {
	my $topLevelDir = shift;
	die "Error lost top level directory" unless defined $topLevelDir;
	my $compressedFastqDir = join('/',$topLevelDir,'fastq');
	return $compressedFastqDir;
}


sub get_sam_mapping_folder {
	my $topLevelDir = shift;
	die "Error lost top level directory" unless defined $topLevelDir;
	my $samDir = join('/',$topLevelDir,'assembly','sam');
	mkdir($samDir, 0777) or die "Error calling mkdir $samDir"
		 unless (-e $samDir);
	return $samDir;
}

sub make_sam_toc {
	my $topLevelDir = shift;
	die "Error lost top level directory" unless defined $topLevelDir;
	my $assemblyDir = shift;
	die "Error lost assemblyDir" unless defined $assemblyDir;
	my $nrOTUs = shift;
	die "Error lost OTU file for bowtie" unless defined $nrOTUs;

	# Where are the gzip'ed FASTQ inputs
	my $compressedFastqDir = get_compressed_fastq_folder($topLevelDir);

	# Where should we find a corresponding SAM file output?
	my $samDir = get_sam_mapping_folder($topLevelDir);

	# make sure there is a SAM file for each input
	# if there isn't one already in place then call bowtie2
	my %input2output;
	my @sams;
        opendir(DIR, $compressedFastqDir) or die "Error opening $compressedFastqDir: $!";
        while(my $file = readdir(DIR)){
		next if ($file =~ /^\./);
		my $basename = basename($file);
		my $input = join('/',$compressedFastqDir,$file);
		my $output = join('/',$samDir,$basename . '.sam');
		$input2output{$input} = $output;
		if (!(-e $output)){
			_bt2_fastq($input,$nrOTUs,$output);
		}else{
			warn "$output already exists so avoiding re-running bt2" if $VERBOSE;
		}
		push(@sams,$output);
        }
        closedir(DIR);

	# make a toc file from all of the sam files
	my $samToc = join('/',$assemblyDir,'sam.toc');
	if (!(-e $samToc)){
		system("$ENV{mPUMA}/bin/create_ace_toc_SAM.pl",
			'-s', 1,	# this is a sorted output, you should be able to disable this if you are chasing speed...
			'-o',	$samToc,
			@sams) == 0 or die "Error creating SAM-based toc file: $!";
	}else{
		warn "$samToc already exists and so avoiding re-parsing all the SAM files" if $VERBOSE;
	}
	return $samToc;
}

sub _bt2_fastq {
	my $queryFile = shift;
	die "Error lost query" unless defined $queryFile;
	my $base = shift;
	die "Error lost base" unless defined $base;
	my $output = shift;
	die "Error lost output" unless defined $output;

	# make sure the base is there
	ensure_bt2_base($base);

	# run bowtie2
	system($bowtie2,
		'--local',
		'-t',
		'-p',	$cpu_cores,
		$base,
		'-U',	$queryFile,
		'-S',	$output,
		'--no-head') == 0 or die "Error calling bowtie2 $queryFile $base $output: $!";
}	

sub ensure_bt2_base {
	my $base = shift;
	die "Error lost bowtie2 base" unless defined $base;
	if (!(
		(-e $base . '.1.bt2')
		and (-e $base . '.2.bt2')
		and (-e $base . '.3.bt2')
		and (-e $base . '.4.bt2')
		and (-e $base . '.rev.1.bt2')
		and (-e $base . '.rev.2.bt2') )){
		# need to build the ref
		system($bowtie2_build,$base,$base)==0 or die "Error calling bowtie2-build $base: $!";
	}
}

sub get_IDs{
        my $file = shift;
        my %IDs;
	if ((!defined $file)or (!(-e $file))){
		warn "WARNING asked to locate IDs from an undefined file. This may be OK if you are skipping chimera checking etc. Otherwise you may be illiciting an unexpected behaviour.";
	}else{	
        	open(FILE,"<$file") or die "Error opening $file for reading : $!";
        	while(<FILE>){
        	        chomp $_;
                	my ($id) = split(/\s+/,$_);
                	$IDs{$id}++;
        	}
        	close FILE or die "Error closing $file: $!";
        }
	return \%IDs;
}

sub make_total_toc {
	my $assemblyDir = shift;
	die "Error lost assemblyDir" unless defined $assemblyDir;
	my $totalReadFile = shift;
	die "Error lost total read set" unless defined $totalReadFile;
	my $aceToc = shift;
	die "Error lost aceToc" unless defined $aceToc;
	my $seqcleanFile = shift;
	die "Error lost seqclean file " unless defined $seqcleanFile;
	my $clusterReport = shift;
	die "Error lost cluster report" unless defined $clusterReport;
	my $chimericIDsFile = shift;
	die "Error lost chimeric IDs file" unless defined $chimericIDsFile;
	my $useChimeraCalls = shift;
        $useChimeraCalls = 1 unless defined $useChimeraCalls;

	# figure out what OTUs were flagged as chimeric
	my $chimericIDs = get_IDs($chimericIDsFile);

	my $totalToc = join('/',$assemblyDir,'total.toc');
	# catch that the file may already exist 
	if (-s $totalToc){
		warn "$totalToc already exists so not regenerating it" if $VERBOSE;
		return $totalToc;
	}
	open(OUTPUT,">$totalToc") or die "Error opening $totalToc for writting: $!";
	print OUTPUT join(",",'ID','nrOTU','Rationale'),"\n";

	# get all of the ReadIDs
	warn "Getting readIDs from $totalReadFile" if ($VERBOSE);
	my $originalReads = get_IDs_from_fasta($totalReadFile);

	# figure out where READs went in the assembly process
	my %readsToInitialOTU;
	warn "Reading initial OTU results from $aceToc" if ($VERBOSE);
	open(FILE,"<$aceToc") or die "Error opening $aceToc for reading: $!";
	while(<FILE>){
		chomp $_;
		my ($otu,$read) = split(/\s+/,$_);
#		next unless ($read eq 'HM2NFJ202JPUMI');
		if (defined $readsToInitialOTU{$read}){
#			warn "$read was placed into more than 1 initial OTU. This is a known feature of gsAssembler in cdna mode. Only tracking the first occurrance";
		}else{
			$readsToInitialOTU{$read} = $otu;
		}
	}
	close FILE or die "Error closing $aceToc: $!";

	# get the OTUs which exist post seqCleaning 
	warn "Reading seqClean results from $seqcleanFile" if $VERBOSE;
	my $postSeqCleanOTUs = get_IDs_from_fasta($seqcleanFile);

	# get the clustering info
	my %nrOTUs;
	my %redundantOTUmapping;
	warn "Reading CD-hit report from $clusterReport" if $VERBOSE;
	open(FILE,"<$clusterReport") or die "Error opening $clusterReport: $!";
	my $header = <FILE>;
	while(<FILE>){
		chomp $_;
		my ($nrOTU,$length,$num,$csv) = split(/\t/,$_);
		$nrOTUs{$nrOTU}++;
		foreach my $rOTU(split(/\,/,$csv)){
			$redundantOTUmapping{$rOTU} = $nrOTU;
		}
	}
	close FILE or die "Error closing $clusterReport: $!";


	warn "Writting output for $totalToc" if $VERBOSE;
	my %finalLocationOfReads;
	foreach my $ID (keys %{$originalReads}){
#		next unless ($ID eq 'HM2NFJ202JPUMI');
#		my $initOTU = $readsToInitialOTU{$ID};
#		my $chimericID = $chimericIDs->{$readsToInitialOTU{$ID}};
#		my $nrOTU = $nrOTUs{$readsToInitialOTU{$ID}};
#		my $rOTU = $redundantOTUmapping{$readsToInitialOTU{$ID}};
#		print <<EOF;
#$ID
#initOTU		$initOTU
#chimericID	$chimericID
#useChimeraCalls	$useChimeraCalls
#nrOTU		$nrOTU
#rOTU		$rOTU
#EOF
#die "HERE";
#	

		my ($finalOTU,$reason) = ('Orphaned','gsAssembler');
		if (defined $readsToInitialOTU{$ID}){
			# read was assembled
			if ($chimericIDs->{$readsToInitialOTU{$ID}}){
				if ($useChimeraCalls){
					# OTU was flagged as possible chimera and we are using those calls
					# OTU was seqcleaned but then trashed because it appeared chimera-like
					$finalOTU	= 'None',
					$reason		= 'Chimeric-OTU',
				}else{
					# OTU was flagged as possible chimera but we are NOT going to make those assertions
					# this would likely be the case when the user has reason to believe that the reference database is lacking for the
					# ecological niche under study
	                                # OTU survived seqCleaning
	                                if ($nrOTUs{$readsToInitialOTU{$ID}}){
	                                        # this OTU was a non-redundant one
	                                        $finalOTU       = $readsToInitialOTU{$ID};
	                                        $reason         = 'CDhit-OTU-was-non-redundant-possible-chimera';
	                                }else{
	                                        $reason         = 'CDhit-OTU-was-redundant-possible-chimera';
	                                        if (defined $redundantOTUmapping{$readsToInitialOTU{$ID}}){
	                                                $finalOTU       = $redundantOTUmapping{$readsToInitialOTU{$ID}};
	                                        }else{
	                                                $finalOTU       = 'Orphaned';
	                                        }
	                                }
				}
			}elsif ($postSeqCleanOTUs->{$readsToInitialOTU{$ID}}){
				# OTU survived seqCleaning
				if ($nrOTUs{$readsToInitialOTU{$ID}}){
					# this OTU was a non-redundant one
					$finalOTU	= $readsToInitialOTU{$ID};
					$reason		= 'CDhit-OTU-was-non-redundant';
				}else{
					$reason		= 'CDhit-OTU-was-redundant';
					if (defined $redundantOTUmapping{$readsToInitialOTU{$ID}}){ 
						$finalOTU	= $redundantOTUmapping{$readsToInitialOTU{$ID}};
					}else{
						$finalOTU	= 'Orphaned';
					}
				}
			}else{
				# did not survive seqClean
				$finalOTU	= 'None';
				$reason		= 'seqClean-trashed';
			}
		}else{
			# read is a singleton
			$finalOTU	= 'Singleton';
			$reason		= 'gsAssembler-left-out';
		}
		print OUTPUT join(",",$ID,$finalOTU,$reason),"\n";
	}
	close OUTPUT or die "Error closing $totalToc: $!";
	return $totalToc;
}


sub get_IDs_from_fasta{
	my $file = shift;
	die "Error lost file" unless defined $file;
	# get all of the ReadIDs
	my %reads;
	my $inio = Bio::SeqIO->new(-file => "<$file", -format => 'fasta') or die "Error opening $file for reading: $!";
	while(my $seq = $inio->next_seq()){
		$reads{$seq->display_id()}++;
	}
	return \%reads;
}



#############################################################################
# post_assembly_cleanup - Runs seqclean and cd_hit_EST clustering on the    #
#                         assembled results.                                #
#############################################################################
# Input Parameters:                                                         #
#                   $assemblyDir - The assembly directory                   #
# Output Parameters:                                                        #
#                   $seqcleanCDHitClustFnaFile - The cluster-cd-hit.fna     #
#						 file.                      #
#                   $seqcleanCDHitClustReportFile - The cluster-cd-hit      #
#						    report file.            #
#                   $seqcleanFile - Non-redundant set of OTU sequences      #
#                                                                           #
#############################################################################

sub post_assembly_cleanup {
	my $assemblyDir = shift;
	die "Error lost assemblyDir" unless defined $assemblyDir;
	my $chimeraCheck = shift;
	die "Error lost whether to use the C3 chimera report or not" unless defined $chimeraCheck;
		my $initialOTUfile = shift;
	die "Error lost initial OTU file" unless defined $initialOTUfile;
	
	# run seqclean to screen isotigs for forward or reverse PCR primers to exclude
	my $isotigsFile = join('/',$assemblyDir, $initialOTUfile);

	# do the seqcleaning
	my $seqcleanFile = seqcleanFile($isotigsFile, $primer_file);

	# insert Chimera checking her ala Chaban's Chimera Checker
	my ($noChimeraFasta,$chimeraIDs) = removePossibleChimeras($seqcleanFile);

	# Cluster the OTU, retreive the longest OTU and make a FASTA of those non-redundant OTU
	my ($seqcleanCDHitClustFnaFile,$seqcleanCDHitClustReportFile);
	if ($chimeraCheck){
		($seqcleanCDHitClustFnaFile,$seqcleanCDHitClustReportFile) = cluster_cd_hit_EST($noChimeraFasta);
	}else{
		# do NOT use the chimera report
		($seqcleanCDHitClustFnaFile,$seqcleanCDHitClustReportFile) = cluster_cd_hit_EST($seqcleanFile);
		$noChimeraFasta = $seqcleanCDHitClustFnaFile;
#		undef $chimeraIDs;
	}

	return ($seqcleanCDHitClustFnaFile, $seqcleanCDHitClustReportFile, $noChimeraFasta, $chimeraIDs);
}
	
sub removePossibleChimeras {
	my $fastaFile = shift;
	die "Error lost fasta file" unless defined $fastaFile;


	# call c2e to get the ends of the sequences
	my $c2eFile = $fastaFile . '.c2e';
	unless (-s $c2eFile){
		system("$ENV{mPUMA}/bin/c2e.pl",
			'-i',	$fastaFile,
			'-o',	$c2eFile,
			'-l', 150) == 0 or die "Error calling c2e.pl on $fastaFile: $!";
	}

	# waterBLAST the ends against the reference DB to see if they are divergent
	my $c2eWateredBlast = $c2eFile . '.wateredBLAST';
	unless (-s $c2eWateredBlast){
		#N.B. this is NOT using the waterBLAST routine within this module because this requires EXACTLY 1 hit and no more for chimera checking
		if ($cpu_cores > 1){
			my $wBcommand = join(' ',$watered_blast,
					'-d', 	$cpnDB_nr_nuc,
					'-p',   'blastn',
					'-h',	0,
					'-v',   1);
                	system("$ENV{mPUMA}/bin/fasta_batcher.pl",
                	        '-i',   $c2eFile,
                	        '-ni',  '-i',
                	        '-o',   $c2eWateredBlast,
                	        '-e',   $wBcommand,
                	        '-c',   $cpu_cores,
                	        '-t',   $cpu_cores,
				'-u',	1) == 0 or die "error calling the wateredblast: $!";
			clean_up_files(abs_path(getcwd()),'^fastaPiece');


		}else{
			system($watered_blast,
				'-i',	$c2eFile,
				'-d',	$cpnDB_nr_nuc,
				'-p',	'blastn',
				'-v',	1,
				'-o',	$c2eWateredBlast) == 0 or die "error waterblast'ing $c2eFile: $!";
		}
	}

	# figure out what is chimeric and what isn't
	my $c3report = $fastaFile . '.c3report';
	unless (-s $c3report){
		system("$ENV{mPUMA}/bin/c3.pl",
			'-i',	$fastaFile,
			'-w',	$c2eWateredBlast,
			'-h',	0,
			'-o',	$c3report) == 0 or die "Error calling c3.pl on $fastaFile: $!";
	}

	# find the IDs of the possible chimeras
	my $chimericIDs = $c3report . '.chimericIDs';
	unless (-s $chimericIDs){
		system("$ENV{mPUMA}/bin/c3extractChimericIDs.pl",
			'-i',	$c3report,
			'-o',	$chimericIDs) == 0 or die "Error finding chimeric IDs from $c3report: $!";
	}

	# make a FASTA which excludes the possible chimeras
	my $noChimeraFile = $fastaFile . '.possibleChimerasRemoved';
	unless(-s $noChimeraFile){
		system("$ENV{mPUMA}/bin/make_fasta_excluding_IDs.pl",
			'-i',	$fastaFile,
			'-e',	$chimericIDs,
			'-o',	$noChimeraFile) == 0 or die "Error creating fasta from $fastaFile while excluding IDs from $chimericIDs: $!";
	}

	return ($noChimeraFile,$chimericIDs);
}

			
sub clean_up_files {
	my $path = shift;
	die "Error lost path" unless defined $path;
	my $pattern = shift;
	die "Error lost pattern" unless defined $pattern;
        my @foundFiles;
        find(sub {
		push(@foundFiles,$File::Find::name) if ($_ =~ /$pattern/);
		},$path);
	
	foreach my $file (@foundFiles){
		unlink $file or die "Error unlinking $file";
	}
}

#############################################################################
# seqcleanFile - Calls seqclean on a Fasta file to remove primers           #
#############################################################################
# Input Parameters:                                                         #
#                   $fasta - The Fasta file to remove primers               #
# Output Parameters:                                                        #
#                   $seqcleanFile - Non-redundant set of OTU sequences      #
#                   $primers - The primer MSF file to lookup primers to     #
#                              remove.                                      #
#                                                                           #
#############################################################################

sub seqcleanFile {
	my $fasta = shift;
	die "Error lost fasta to seqclean" unless defined $fasta;
	
	my $primers = shift;
	$primers = $primer_file unless defined $primers;
	my $seqcleanFile = $fasta . '.seqclean';

	my $seqcleanCmd = "$seqclean $fasta -c $cpu_seqclean -v $primers -o " . $seqcleanFile;
	unless(-s $seqcleanFile){
		# format the msf database file into .nin .nsq .nhr files.
		formatdbNuc($primers);

		# Force seqclean to work where the fasta is so that we don't pollute the cwd with all the cleaning files

		# figure out where we currently are
		my $prev_dir = abs_path(getcwd());

		my $dir = dirname($fasta);
		chdir($dir) or die "Error changing directory to $dir: $!";

		# seqclean the isotigs of the qPCR primers.
		system($seqclean, $fasta, '-c', $cpu_seqclean, '-v', $primers, '-o', $seqcleanFile) == 0 or die "Error calling $seqcleanCmd: $?";	

		# this would be a great spot for some parallelism but seqclean is hardcoded to do a couple of things which REALLY make it tricky. 
		# it HAS to have the input Fasta as ARGV[0]
		# each seqclean process creates a cleaning_CPUNUM directory in cwd to break the chunks into 
		# this means that you have collisions is 2 processes are running in the same directory...
		# Arrrrgh!

		chdir($prev_dir) or die "Error changing directory back to $prev_dir: $!";
	}

	return $seqcleanFile;
}

sub formatdbNuc{
	my $fastadb = shift;
	die "Error lost fasta to seqclean" unless defined $fastadb;
	# format the database file into .nin .nsq .nhr files.
	my ($dbFileNIN, $dbFileNSQ, $dbFileNHR);
	$dbFileNIN = $fastadb . ".nin";
	$dbFileNSQ = $fastadb . ".nsq";
	$dbFileNHR = $fastadb . ".nhr";
	unless(-s $dbFileNIN and -s $dbFileNSQ and -s $dbFileNHR){
		warn "Calling formatdb for $fastadb....\n\n$formatdb -p F -i $fastadb" if $VERBOSE;
		system($formatdb, 
			'-p', 'F', 
			'-i', $fastadb
		) == 0 or die "Error calling $formatdb -p F -i $fastadb: $?";
	}
}

sub formatdbPep{
	my $fastadb = shift;
	die "Error lost fasta to seqclean" unless defined $fastadb;
	# format the database file into .pin .psq .phr files.
	my ($dbFilePIN, $dbFilePSQ, $dbFilePHR);
	$dbFilePIN = $fastadb . ".pin";
	$dbFilePSQ = $fastadb . ".psq";
	$dbFilePHR = $fastadb . ".phr";
	unless(-s $dbFilePIN and -s $dbFilePSQ and -s $dbFilePHR){
		warn "Calling formatdb for $fastadb....\n\n$formatdb -p T -i $fastadb" if $VERBOSE;
		system($formatdb, 
			'-p', 'T', 
			'-i', $fastadb
		) == 0 or die "Error calling $formatdb -p T -i $fastadb: $?";
	}
}

#############################################################################
# cluster_cd_hit_EST - Performs the cd_hit clustering on the fasta file to  #
# obtain the non-redundant OTU seqeunces file                               #
#############################################################################
# Input Parameters:                                                         #
#                   $fasta - The Fasta file to remove primers               #
# Output Parameters:                                                        #
#                   $seqcleanCDHitClustFnaFile - The cluster-cd-hit.fna     #
#						 file.                      #
#                   $seqcleanCDHitClustReportFile - The cluster-cd-hit      #
#						    report file.            #
#                                                                           #
#############################################################################

sub cluster_cd_hit_EST {
	my $fasta = shift;
	die "Error lost fasta" unless defined $fasta;

	# Clustering via CD-hit
	my $CDHitFile = $fasta . '.cd-hit';
	unless(-s $CDHitFile){
		my $cdhitCmd = "$cdHitEST -c 1 -G 1 -n 10 -M 0 -i $fasta -o $CDHitFile";
		warn "$cdhitCmd\n" if $VERBOSE;
		system($cdHitEST, 
			'-c',	1,
			'-G',	1,
			'-n',	10,
			'-M',	0,
			'-d',	0,	# include really long seqIDs 
			'-i',	$fasta, 
			'-o',	$CDHitFile
		) == 0 or die "Error calling $cdhitCmd: $?";
	}

	# Create a summary report from the CD-hit step
	my $CDHitClustFile = $CDHitFile . '.clstr';
	my $CDHitClustReportFile = $CDHitClustFile . '.report';
	unless(-s $CDHitClustReportFile){
		warn $parseCDHitEst . " -c $CDHitClustFile -o $CDHitClustReportFile\n\n" if($VERBOSE);
		system($parseCDHitEst, '-c', $CDHitClustFile, '-o', $CDHitClustReportFile) == 0 or die "Error calling $parseCDHitEst: $?";
	}

	# make sure that the ed results are indexed with CDBFASTA
	my $fnaIdxFile = $fasta . '.cidx';
	unless(-s $fnaIdxFile){
		# Create a CDBFasta index of the ed isotigs
		my $cdbfastaCmd = "$cdbfasta $fasta";
		warn "$cdbfastaCmd\n" if($VERBOSE);
		system($cdbfasta, $fasta) == 0 or die "Error calling $cdbfastaCmd: $?";
	}

	# Create a fasta file of the ed isotigs which were clustered down into a non-redundant set.
	my $CDHitClustFnaFile = $CDHitClustFile . '.fna';
	unless (-s $CDHitClustFnaFile){
		open INFILE, "<$CDHitClustReportFile" or die $!;
		my @isotigs = ();
		my $header = <INFILE>;
		while(<INFILE>){
			chomp $_;
			my ($id) = split(/\t/,$_);
			push(@isotigs, $id);
		}
		close INFILE or die "Error closing $CDHitClustReportFile: $!";

		@isotigs = sort(@isotigs);

		open OUTFILE, ">", $CDHitClustFnaFile or die $!;
		foreach my $id (@isotigs){
			my $cdbyankCmd = "$cdbyank $fnaIdxFile -a '" . $id . "'";
			warn "$cdbyankCmd\n" if $VERBOSE;
		
			local (*CDBYANK_OUT, *CDBYANK_IN);
			my  $pid = open2(\*CDBYANK_OUT,\*CDBYANK_IN, $cdbyankCmd) or die "Error calling open2: $!";
			close CDBYANK_IN or die "Error closing STDIN to cdbyank process: $!";	
			
			while (<CDBYANK_OUT>){ 
				print OUTFILE $_;
			}
			
			close CDBYANK_OUT or die "Error closing STDOUT from cdbyank process: $!";
			wait;
		}
		close OUTFILE or die $!;
	}
	return ($CDHitClustFnaFile,$CDHitClustReportFile);
}

#############################################################################
# translationOfStranded -  Using the BLASTX results solve for the           #
# translation of each isotigs                                               #
#############################################################################
# Input Parameters:                                                         #
#                   $FnaFile - The path to the fna file.                    #
#                   $FnaBlastxFile - The path to the                        #
#                                     fna.blastx file.                      #
# Output Parameters:                                                        #
#                   $FnaAaFile - The path to the fna.aa file.               #
#############################################################################
# 2.	Using the BLASTX results solve for the translation of each isotigs
# a.	$mPUMA/bin/translate_stranded.pl -i 454Isotigs.fna.seqclean -o 454Isotigs.fna.seqclean.aa -b 454Isotigs.fna.seqclean.blastx

sub translationOfStranded{
	my $FnaFile = shift;
	die "Error lost fna file" unless defined $FnaFile;
	my $FnaBlastxFile = $FnaFile . '.blastx';
	
	my $orientedNucFile = $FnaFile . '.oriented';
#	my $FnaAaFile = $FnaFile . '.aa';
	my $FnaAaFile = $orientedNucFile . '.aa';
	# convert the files
	unless(-s $FnaAaFile){
		warn "\nGenerating amino acid fna file....\n" if($VERBOSE);
		my $transCmd  = "$translateStranded -i $FnaFile -o $FnaAaFile -b $FnaBlastxFile\n\n";
		warn $transCmd if $VERBOSE;
		system($translateStranded, # this needs to output an oriented FASTA of the nucleotides too...
			'-i', $FnaFile,
			'-o', $FnaAaFile,
			'-n', $orientedNucFile,
			'-b', $FnaBlastxFile
		) == 0 or die "Error calling $translateStranded: $?";
	}
	return ($orientedNucFile,$FnaAaFile);
}

# takes a fasta file as an argument and does the cd_hit clustering and retuns the path to the fasta of the non-redundant OTU seqeunces
#############################################################################
# cluster_cd_hit - Performs the cd_hit clustering on the aa file to         #
# obtain the non-redundant OTU protein seqeunces file                       #
#############################################################################
# Input Parameters:                                                         #
#                   $FnaAaFile - The path to the fna.aa file.               #
# Output Parameters:                                                        #
#                   $FnaAaCDHitClustFile - The aa.cluster-cd-hit            #
#				           file.                            #
#                                                                           #
#############################################################################

sub cluster_cd_hit {
	my $AaFnaFile  = shift;
	die "Error lost fasta" unless defined $AaFnaFile ;
# a.	/usr/local/cd-hit/cd-hit -c 1 -G 1 -n 5 -M 0 -i 454Isotigs.fna.seqclean.aa -o 454Isotigs.fna.seqclean.aa.cd-hit
	# Clustering via CD-hit
	my $AaCDHitFnaFile = $AaFnaFile . '.cd-hit';
	unless(-s $AaCDHitFnaFile){
		my $cdhitCmd = "$cdHit -c 1 -G 1 -n 5 -M 0 -i $AaFnaFile  -o $AaCDHitFnaFile";
		warn "$cdhitCmd\n" if $VERBOSE;
		system($cdHit, 
			'-c',	1,
			'-G',	1,
			'-n',	5,
			'-M',	0,
			'-i',	$AaFnaFile , 
			'-o',	$AaCDHitFnaFile
		) == 0 or die "Error calling $cdhitCmd: $?";
	}

	# Create a summary report from the CD-hit step
	my $AaCDHitClustFile = $AaCDHitFnaFile . '.clstr';
	my $AaCDHitClustReportFile = $AaCDHitClustFile . '.report';
	unless(-s $AaCDHitClustReportFile){
		warn $parseCDHitEst . " -c $AaCDHitClustFile -o $AaCDHitClustReportFile\n\n" if($VERBOSE);
		system($parseCDHitEst, 
			'-c', $AaCDHitClustFile, 
			'-o', $AaCDHitClustReportFile
		) == 0 or die "Error calling $parseCDHitEst: $?";
	}
	
		# make sure that the ed results are indexed with CDBFASTA
	my $fnaIdxFile = $AaFnaFile  . '.cidx';
	unless(-s $fnaIdxFile){
		# Create a CDBFasta index of the ed isotigs
		my $cdbfastaCmd = "$cdbfasta $AaFnaFile ";
		warn "$cdbfastaCmd\n" if($VERBOSE);
		system($cdbfasta, $AaFnaFile ) == 0 or die "Error calling $cdbfastaCmd: $?";
	}

	# Create a fasta file of the ed isotigs which were clustered down into a non-redundant set.
	my $AaCDHitClustFnaFile = $AaCDHitClustFile . '.fna';
	unless (-s $AaCDHitClustFnaFile){
		open INFILE, "<$AaCDHitClustReportFile" or die $!;
		my @isotigs = ();
		my $header = <INFILE>;
		while(<INFILE>){
			chomp $_;
			my ($id) = split(/\t/,$_);
			push(@isotigs, $id);
		}
		close INFILE or die "Error closing $AaCDHitClustReportFile: $!";

		@isotigs = sort(@isotigs);

		open OUTFILE, ">", $AaCDHitClustFnaFile or die $!;
		foreach my $id (@isotigs){
			my $cdbyankCmd = "$cdbyank $fnaIdxFile -a '" . $id . "'";
			warn "$cdbyankCmd\n" if $VERBOSE;
		
			local (*CDBYANK_OUT, *CDBYANK_IN);
			my  $pid = open2(\*CDBYANK_OUT,\*CDBYANK_IN, $cdbyankCmd) or die "Error calling open2: $!";
			close CDBYANK_IN or die "Error closing STDIN to cdbyank process: $!";	
			
			while (<CDBYANK_OUT>){ 
				print OUTFILE $_;
			}
			
			close CDBYANK_OUT or die "Error closing STDOUT from cdbyank process: $!";
			wait;
		}
		close OUTFILE or die $!;
	}
	return ($AaCDHitClustFnaFile, $AaCDHitClustReportFile);
}

sub t_coffee_fasta {
	my $fastaFile = shift;
	die "Error lost file" unless defined $fastaFile;

	my $output = $fastaFile . '.aln';
	unless (-s $output){

                # figure out where we currently are
                my $prev_dir = abs_path(getcwd());

                my $dir = dirname($fastaFile);
                chdir($dir) or die "Error changing directory to $dir: $!";

		system($tCoffee,
			$fastaFile,
			'-outfile',	$output,
			'-mode',	'quickaln',
			'-output',	'fasta_aln') == 0 or die "Error calling $tCoffee: $!";

                chdir($prev_dir) or die "Error changing directory back to $prev_dir: $!";
	}
	return $output;
}

sub fastTree_nuc {
	my $alnFile = shift;
	die "Error lost alnFile" unless defined $alnFile;
	my $useGTR = shift;
	$useGTR = 0 unless defined $useGTR;
	my @params = ('-nt');
	push(@params,'-gtr') if $useGTR;
	return fastTree_aln($alnFile,\@params);
}

sub fastTree_aln {
	my $alnFile = shift;
	die "Error lost alnFile" unless defined $alnFile;
	my $params = shift;

	# N.B. this could be modified to change the BLOSUM matrix default to use the cpn60 UT specific matrix... 
	# need to find a way to confirm how to evaluate the two trees as to which is better...

	my $output = $alnFile . '.FastTree';
	unless (-s $output){
		open(OUT,">$output") or die "Error opening $output for writting: $!";
		warn "Calling $fastTree ",join(' ',@{$params})," $alnFile" if ($VERBOSE and defined $params);
		
		my $pid = open2(\*CHLDOUT,\*CHLDIN,$fastTree,'-fastest',@{$params},$alnFile) 
			or die "Error calling open2 on $fastTree: $!";
		close CHLDIN or die "Error closing STDIN to child process: $!";
		print OUT while(<CHLDOUT>);
		close CHLDOUT or die "Error closing STDOUT from child process: $!";
		close OUT or die "Error closing $output: $!";
		wait;
	}
	return $output;
}


#############################################################################
# waterBlastIT - Runs a waterBlast on a query and target DNA or RNA         #
# ( BLASTN ) sequences to obtain a similarity comparison of the sequencies. #
#############################################################################
# Input Parameters:                                                         #
#                   $wbQuery - The 454Isotigs.fna query file                #
#                   $wbTarget - The cpndb_nr database target file           #
#		    blast_prog - must be either blastp or blastn.           #
# Output Parameters:                                                        #
#                   $wbSymlinkFile - waterBlast format output symbolic      #
#                                    link file                              #
#                   $status - Status of the waterBlast check                #
#                                                                           #
#############################################################################

# Changes ~ July 7th 
# simple script to return what the actual target file should be ...
# 
sub nameWaterBlastResultFile {
	my $wbQuery = shift;
	die "Error lost the 454Isotigs.fna query file" unless defined $wbQuery;
	my $wbDB = shift;
	my ($dbName,$dbPath) = fileparse($wbDB);
	return $wbQuery . '.wateredBLAST.' . $dbName;
}

sub waterBlastIT{
	my $wbQuery = shift;
	die "Error lost the query file" unless defined $wbQuery;
	my $wbDB = shift;
	die "Error lost the database target file" unless defined $wbDB;
	my $blast_prog = shift;
        die "Error lost the blast program ( blastp or blastn )" unless defined $blast_prog;

	warn "\nSequence DB comparison of the assembled sequences....\nProcessing Job: waterBlast on $wbQuery and $wbDB \n\n" if ($VERBOSE);
	my $wbOutputFile = nameWaterBlastResultFile($wbQuery,$wbDB);

	if (-s $wbOutputFile) {
		# Checks the waterBlast output file for errors.
		my $status = WaterBlastCheck($wbOutputFile);
		return ($wbOutputFile, $status);
 	} 

	# Note how this command an arguments are being passed to the open2 call. This is much safer. 
	# wateredBlast I/O command
	# Don't need open2 for this 

	if ($cpu_cores > 1){
		my $wbcommand = join(' ',$watered_blast,
				'-v',   $num_descs,
				'-p',   $blast_prog,
				'-d',   $wbDB,
				'-h',	0);
               	system("$ENV{mPUMA}/bin/fasta_batcher.pl",
               	        '-i',   $wbQuery,
               	        '-ni',  '-i',
               	        '-o',   $wbOutputFile,
               	        '-e',   $wbcommand,
               	        '-c',   $cpu_cores,
               	        '-t',   $cpu_cores
		) == 0 or die "error calling the wateredblast: $!";
		clean_up_files(abs_path(getcwd()),'^fastaPiece');

	}else{
		system($watered_blast,
			'-v',	$num_descs,
			'-p',	$blast_prog,
			'-i',	$wbQuery,
			'-d',	$wbDB,
			'-o',	$wbOutputFile,
			'-h',	0
		) == 0 or die "Error calling watered_blast: $!";


	}

	# Checks the waterBlast output file for errors.
	my $status = WaterBlastCheck($wbOutputFile);
	return ($wbOutputFile, $status);
}

#############################################################################
# WaterBlastCheck - Checks the waterBlast output file for errors.           #
# If success returns a sucessful message                                    #
# Otherwise kills program and prints the file/line number of the error.     #
############################################################################
# Input Parameters:                                                         #
#		    $wbOutputFileName - waterBlast format output file	    #
#				        to check		            #
# Output Parameters:                                                        #
#                   $status - Status of the waterBlast check                #
#                                                                           #
#############################################################################

sub WaterBlastCheck{

	my ($failure, $error) = "";
	my $status = -1;
	my $i = 1;

	my $wbOutputFileName = shift;
	die "Error lost the waterBlast format output file to check" unless defined $wbOutputFileName;

	open(INFILE, "<$wbOutputFileName") or die "Error opening $wbOutputFileName for reading: $!"; #open for read
	# N.B. this assumes no header!!!

	while(<INFILE>){
		my @line = ();
		@line = split("\t", $_);
		if(scalar(@line) == 5 or scalar(@line) == 6){
			if(($line[2] =~ m/\d+.\d+|0|N\/A/) and ($line[3] =~ m/\d+|N\/A/) and ($line[4] =~ m/.+|N\/A/) ){
				$status = 1;
			}
			else{
				$status = 0;
				if ($status == 0){
					$failure = "failure";
				}
				$error = "Error in file: " . $wbOutputFileName . " line# " . $i . " status: " . $failure . "\n" . $_;
				$status = $error;
				return $status;
			}
		}
		else{
			$status = 0;
			if ($status == 0){
				$failure = "failure";
				$error = "Error in file: " . $wbOutputFileName . " line# " . $i . " status: " . $failure . "\n" . $_;
				$status = $error;
				return $status;
			}
		}
		$i++;
	}

	close INFILE or die $!;
	return $status;
}

#############################################################################
# solveOrientation - Generates the stranded fasta file and returns the path #
#                    to that file                                           #
#############################################################################
# Input Parameters:                                                         #
#                   $fasta - The Fasta file to remove primers               #
#                   $wbOutput - waterBlast format output file               #
# Output Parameters:                                                        #
#		    $stranded - stranded fasta format output file	    #
#                                                                           #
#############################################################################

sub solveOrientation {
	my $fasta = shift;
	die "Error lost fasta" unless defined $fasta;

	my $wbOutput = shift;
	die "Error lost wateredBLAST output" unless defined $wbOutput;

	my $stranded = $fasta . '.stranded';

	# Create a properly oriented  / stranded fasta of the cpn60 universal targets based on the wateredBLAST report
	unless(-s $stranded){
		my $fixStrandCmd = $fixStrand . " -i $fasta -w $wbOutput -o $stranded";
		warn "\n$fixStrandCmd\n" if $VERBOSE;
		system($fixStrand, '-i', $fasta, '-w', $wbOutput, '-o',  $stranded) == 0 or die "Error calling $fixStrandCmd: $?";
	}
	return $stranded;
}

#############################################################################
# sffinfoFna - Creates a total read set file "total.fna" from the SFF       #
# files.                                                                    #
#############################################################################
# Input Parameters:                                                         #
# 		    $dir - The directory to store the "total.fna" file	    #
#                                                                           #
# Output Parameters:                                                        #
#		    $totalReadSet - "total.fna" format output file	    #
#                   $status - Status of the sffinfoFna check                #
#                                                                           #
#############################################################################

sub sffinfoFna{
	my $dir = shift;
	die "error lost dir" unless defined $dir;
	my $totalReadSet = join('/',$dir,'total.fna');
	my $status = 0;
	#@_ is the list of SFFs
	warn "\nGenerating total read set file \"$totalReadSet\"....\n\n" if ($VERBOSE);
	if (-s $totalReadSet) {

 		warn "$totalReadSet exists!....\n" if($VERBOSE);

		# Checks the format of the total readset file "total.fna" output file for errors.
		$status = sffinfoFnaCheck($totalReadSet);
		return ($totalReadSet, $status);
 	} 

	# make the FNA file
	$totalReadSet = sffsToFna($totalReadSet,@_);	
	# Checks the format of the total readset file "total.fna" output file for errors.
	$status = sffinfoFnaCheck($totalReadSet);
	
	return ($totalReadSet, $status);
	
}

#############################################################################
# sffsToFna - Creates a total read set file "total.fna" from the SFF files. # 
# Calls sffinfo on all of the SFF files and writes the result to $output    #
#############################################################################
# Input Parameters:                                                         #
#		    $output - Name and filepath of the "total.fna" output   #
#                             file.                                         #
#                                                                           #
# Output Parameters:                                                        #
#		    $output - Name of the "total.fna" output file.          #
#                                                                           #
#############################################################################

sub sffsToFna {
	my $output = shift;
	die "Error lost file to output to" unless defined $output;

	local (*SFFINFO_OUT, *SFFINFO_IN);
	open (OUTFILE, ">$output") or die "Error opening $output for writting: $!";
	foreach my $listingSFFs (@_){
		my $pid = open2(\*SFFINFO_OUT,\*SFFINFO_IN,$sffinfo,'-s',$listingSFFs) or die "Error calling open2: $!";
		close SFFINFO_IN or die "Error closing STDIN to sffinfo_fna process: $!";	#have to close Writer before read		
		while (<SFFINFO_OUT>){ 
			print OUTFILE $_;
		}
		close SFFINFO_OUT or die "Error closing STDOUT from sffinfo_fna process: $!";
		wait;
	}
	close OUTFILE or die $!;

	return $output;
}

#############################################################################
# sffinfoFnaCheck - Checks the "total.fna" output file for errors.	    #
# If success returns successful message                                     #
# Otherwise kills program and prints the file/line number of the error.     #
#############################################################################
# Input Parameters:                                                         #
#		    $file to check				            #
#                                                                           #
# Output Parameters:                                                        #
#                   $status - Status of the sffinfoFna check                #
#                                                                           #
#############################################################################

sub sffinfoFnaCheck{

	my $file = shift;
	die "Error lost file" unless defined $file;

	my ($failure, $error) = "";
	my $status = -1;
	open (INFILE, "<$file") or die "Error opening $file for reading: $!";

	my $i = 1;
 	while (<INFILE>){
		chomp $_;
		next if ($_ =~ /^>/);
		if($_ =~ /^[ACGTN]+$/i){
			$status = 1;
		}
		else{
			$status = 0;
			if ($status == 0){
				$failure = "format failure";
				$error = "Error in file: \"" . $file . "\" line# " . $i . " status: " . $failure . "\n" . $_;
				$status = $error;
				return $status;
			}
		}
		$i++;
	}
	close INFILE or die $!;
	return $status;
}

#############################################################################
# calcSpSnMetrics - Computes the Specificity and Sensitivity on a query     #
# and target sequence to obtain metrics validating the waterBlast hits.     #
#############################################################################
# Input Parameters:                                                         #
#		    $totalReadSet - The total.fna query file		    #
# 		    $contigFileName - The 454Isotigs.fna query file	    #
# 		    $assemblyDir - The assembly directory		    #
#                                                                           #
# Output Parameters:                                                        #
#		    $snspInputFile - Input for the calSpSnMetricsCheck sub  #
# 				     routine 		                    #
#				     (unaltered "Specificity_Sensitivity"   #
# 				     file)		                    #
#                                                                           #
#		    $snspOutputFile - Input for the calSpSnMetricsCheck     #
# 				      sub routine 		            #
# 				      (empty "Specificity_Sensitivity"      #
# 				      file for extracting all useful	    #
# 				      information separately.)		    #
#                   $status - Status of the calcSpSnMetrics check           #
#                                                                           #
#############################################################################

sub calcSpSnMetrics{

	my $assemblyDir = shift;
	die "Error lost assembly dir" unless defined $assemblyDir;
	my $totalReadSet = shift;
	die "Error lost totalReadSet" unless defined $totalReadSet;
	my $contigFileName = shift;
	die "Error lost totalReadSet" unless defined $contigFileName;
	my $totalReadWB = shift;
	die "Error lost watered blast for total " unless defined $totalReadWB;
	my $contigWB = shift;
	die "Error lost contig watered blast" unless defined $contigWB;
	warn "\nGenerating validation metrics of the assemblies....\n\n" if ($VERBOSE);
	my $snspInputFile = join('/',$assemblyDir,'Specificity_Sensitivity');
	my $snspOutputFile = join('/',$assemblyDir,'calSpSnMetrics');
	$contigFileName =~ s/$assemblyDir//g;
	warn "calculate_Sn_Sp_metrics on the assembly $assemblyDir within input from $totalReadSet and $contigFileName using reference DB $cpnDB_nr_nuc\n\n" if ($VERBOSE);

	if (-s $snspOutputFile) {
 		warn "$snspInputFile exists!....\n" if ($VERBOSE);
		# Checks the "Specificity_Sensitivity" output file for errors.
		my $status = calcSpSnMetricsCheck($snspInputFile, $snspOutputFile);
		return ($snspInputFile, $snspOutputFile, $status);
 	} 
	my $calSpSnMetricsCmd =  "$calculate_sp_sn_metrics -a $assemblyDir -i $totalReadSet -v $num_descs -c $contigFileName -d $cpnDB_nr_nuc -verbose 1 -s 1 -o $snspInputFile -iw $totalReadWB -cw $contigWB\n\n";
	warn $calSpSnMetricsCmd if ($VERBOSE);
	# calSpSnMetrics I/O command
	system($calculate_sp_sn_metrics,
		'-a',	$assemblyDir,
		'-i',	$totalReadSet,
		'-v',	$num_descs,
		'-d',	$cpnDB_nr_nuc,
		'-c',	$contigFileName,
		'-verbose',	1,
		'-s',	1,				# the watered blast searches are already done and in place
		'-o',	$snspInputFile,			# Added an output parameter to the script - this needs to be tested to make sure it works correctly
		'-iw',	$totalReadWB,
		'-cw',	$contigWB
	) == 0 or die "Error calling $calculate_sp_sn_metrics: $!";
	
	# Checks the "Specificity_Sensitivity" output file for errors.
	my $status = calcSpSnMetricsCheck($snspInputFile, $snspOutputFile);
	return ($snspInputFile, $snspOutputFile, $status);

}

#############################################################################
# calcSpSnMetricsCheck - Checks the "Specificity_Sensitivity" output file   #
# for errors.                                                               #
# If success returns successful message                                     #
# Otherwise kills program and prints the file/line number of the error.     #
#############################################################################
# Input Parameters:                                                         #
#		    $snspInputFile - Input for the calSpSnMetricsCheck sub  #
# 				     routine. 	    			    #
#				     (unaltered "Specificity_Sensitivity"   #
#				     file)		    		    #
#		    $snspOutputFile - Input for the calSpSnMetricsCheck     #
# 				      sub routine.		    	    #
#				      (empty "Specificity_Sensitivity"	    #
# 				      file for extracting all useful	    #
# 				      information separately will	    #
# 				      exclude all N/A lines). 		    #
#                                                                           #
# Output Parameters:                                                        #
#                   $status - Status of the calcSpSnMetrics check           #
#                                                                           #
#############################################################################

sub calcSpSnMetricsCheck{	

	my ($failure, $error) = "";
	my $status = -1;

	my $snspInputFile = shift;
        die "Error lost the unaltered \"Specificity_Sensitivity\" file" unless defined $snspInputFile;
	my $snspOutputFile = shift;
        die "Error lost the empty \"calcSpSnMetrics\" file" unless defined $snspOutputFile;

	open(INFILE, "<", "$snspInputFile") or die $!; #open for read
	open(OUTFILE, ">", "$snspOutputFile") or die $!; #open for write
	my $i = 1;
	while(<INFILE>){
		if($i ne 1){
			if($_ =~ m/.+\t\d+.\d+|0|N\/A\t\d+.\d+|0|N\/A\t\d+.\d+|0|N\/A\t\d+\t.+\n/){
				$status = 1;
				if($_ !~ m/.+\tN\/A\tN\/A\tN\/A\t\d+\tUNKNOWN\n/){
					print OUTFILE $_;
				}
			}
			else{
				$status = 0;
				if ($status == 0){
					$failure = "format failure";
					$error = "Error in file: " . $snspInputFile . " line# " . $i . " status: " . $failure . "\n" . $_;
					$status = $error;
					return $status;
				}
			}
			
		}	
		elsif($i eq 1){
			warn "\nChecking the calSpSnMetrics output file " . $snspInputFile . " for errors...\n" if ($VERBOSE);
			print OUTFILE $_;
		}
		$i++;
	}
	close INFILE or die $!;
	close OUTFILE or die $!;
	return $status;
}

#############################################################################
# createLibraries - Creates a libraries dir and then a sub folder for each  #
#                   library.                                                #
#############################################################################
# Input Parameters:                                                         #
#                   $dir - The working directory                            #
#                   $lib_info  - Hash reference pointing to the SFF         #
#                                libraries.                                 #
# Output Parameters:                                                        #
#                   $libraries_dir - The libraries directory                #
#                                                                           #
#############################################################################

sub createLibraries{
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	my $lib_info = shift;
	die "Error lost library info" unless defined $lib_info;

	# what is the sub-folder which will hold all of the libraries
	my $libraries_dir = join('/',$dir,'libraries');
	unless(-d $libraries_dir){
		mkdir($libraries_dir, 0777) or die "Can't make directory $libraries_dir: $!";
	}

	foreach my $lib (keys %{$lib_info}){
		my $d = join('/',$libraries_dir,$lib);
		next if (-d $d);
		mkdir($d, 0777) or die "Can't make directory $d: $!";
	}		
	return $libraries_dir;
}

#############################################################################
# generateReadsIds - Creates the "reads.ids" files for each subsequent SFF  #
# directory in the library directory.                                       #
#############################################################################
# Input Parameters:                                                         #
#                   $libDir - The library directory                         #
#                   $sffFiles - The SFF filename path                       #
# Output Parameters:                                                        #
#                   $numReads - The number of reads in the "reads.ids"      #
#                               file                                        #
#                                                                           #
#############################################################################

sub generateReadIds{

	my $libDir = shift;
	die "Error lost library dir" unless defined $libDir;

	my $sffFiles = shift;
	die "Error lost sff file list" unless defined $sffFiles;

	my $readsIdsOutputFile = join('/',$libDir,$read_ids_file);	# $read_ids_file is defined in this module

	my $numReads = 0;

	if (-s $readsIdsOutputFile){
		warn "$readsIdsOutputFile file exists! refusing to write to this file in case this would lead to duplicating a set of IDs....\n" if ($VERBOSE);
		open(READS,"<$readsIdsOutputFile") or die "Error opening $readsIdsOutputFile for reading: $!";
		while(<READS>){
			chomp $_;
			$numReads++;
		}
		close READS or die "Error closing $readsIdsOutputFile: $!";
	}else{
		open(OUTFILE,">$readsIdsOutputFile") or die "Error opening $readsIdsOutputFile for writting: $!";
		foreach my $sff (@{$sffFiles}){
			local (*SFFINFO_OUT, *SFFINFO_IN);
			my $pid = open2(\*SFFINFO_OUT,\*SFFINFO_IN, $sffinfo, '-a', $sff) or die "Error calling open2: $!";
			close SFFINFO_IN or die "Error closing STDIN to sffinfo_fna process: $!";	#have to close Writer before read
			while (<SFFINFO_OUT>){ 
				print OUTFILE $_;
				$numReads++;
			}
			close SFFINFO_OUT or die "Error closing STDOUT from sffinfo_fna process: $!";
			wait;
		}
		close OUTFILE or die "Error closing $readsIdsOutputFile: $!";
		# There is no check on the read ids worht doing unless you are checking the number of lines and comparing to $numReads from here.
	}

	return $numReads;
}


sub generateReadIdsFastq{
	my $topLevelDir = shift;
	die "Error lost top level directory" unless defined $topLevelDir;

	my $libDir = shift;
	die "Error lost library dir" unless defined $libDir;

	my $ngsFiles = shift;
	die "Error lost file list" unless defined $ngsFiles;

        my $compressedFastqDir = get_compressed_fastq_folder($topLevelDir);

	my $readsIdsOutputFile = join('/',$libDir,$read_ids_file);	# $read_ids_file is defined in this module

	my $numReads = 0;

	if (-s $readsIdsOutputFile){
		warn "$readsIdsOutputFile file exists! refusing to write to this file in case this would lead to duplicating a set of IDs....\n" if ($VERBOSE);
		open(READS,"<$readsIdsOutputFile") or die "Error opening $readsIdsOutputFile for reading: $!";
		while(<READS>){
			chomp $_;
			$numReads++;
		}
		close READS or die "Error closing $readsIdsOutputFile: $!";
	}else{
		open(OUTFILE,">$readsIdsOutputFile") or die "Error opening $readsIdsOutputFile for writting: $!";
		foreach my $originalFile (@{$ngsFiles}){
			my $basename = basename($originalFile);
			my $target;
			if ($basename =~ /.*\.gz$/){
				# original file was compressed
				$target = join('/',$compressedFastqDir,$basename);
			}else{
				# original file was not compressed
				$target = join('/',$compressedFastqDir,$basename.'.gz');
			}
			if (-e $target){
				# read the file to find the # of reads and output the IDs
				my $in = IO::Uncompress::Gunzip->new($target) or die "Error opening $target for reading: $!";	
				my $inio = Bio::SeqIO->new(-Fh => $in, -format => 'fastq-illumina');
				while(my $seq = $inio->next_seq()){
					$numReads++;
					print OUTFILE $seq->display_id(),"\n";
				}
				undef $inio;
				undef $in;
			}else{
				die "Error cannot figure out the fastq file for $originalFile => $basename in $compressedFastqDir";
			}
		
		}
		close OUTFILE or die "Error closing $readsIdsOutputFile: $!";
		# There is no check on the read ids worht doing unless you are checking the number of lines and comparing to $numReads from here.
	}

	return $numReads;
}


sub check_for_ace_compression{
	my $assemblyDir = shift;
	die "Error lost assemblyDir" unless defined $assemblyDir;

	# Find any ace files in the assembly dir
	my @aceFiles;
	find(sub {
			if ($_ =~ /\.ace$/){
				push(@aceFiles,$File::Find::name);
			}
		},$assemblyDir);

	foreach my $file (@aceFiles){
		warn "Trying to compress $file" if $VERBOSE;
		system('gzip',$file) == 0 or die "Error compressing $file: $!";
	}
}

sub run_classifier_on_fasta {
	my $input = shift;
	die "Error lost input" unless defined $input;
	my $output = $input . '.classified';

	if (-s $output){
		warn "$output already exists so skipping" if $VERBOSE;
		return $output;
	}else{
		# cpn60
		my $javaClassifierCommand = "java -Xmx1024m -jar $classifier -t $ENV{mPUMA}/classifier/rRNAClassifier.properties";
		# 16S rRNA
#		my $javaClassifierCommand = "java -Xmx1024m -jar $classifier";
		system("$ENV{mPUMA}/bin/fasta_batcher.pl",
			'-i',	$input,
			'-ni',	'-q',
			'-o',	$output,
			'-e',	$javaClassifierCommand,
			'-c',	$cpu_cores,
			'-t',	$cpu_cores) == 0 or die "Error calling the classifier"; 
		return $output;
	}
}

sub subsample_library {
	my $libDir = shift;
	die "Error lost library dir" unless defined $libDir;
	my $number = shift;
	die "Error lost subsampling depth" unless defined $number;

        my $src = join('/',$libDir,$read_ids_file);
        my $dest = join('/',$libDir,$read_sub_file);
	if (-e $dest){
		warn "Not downsampling as $dest already exists" if $VERBOSE;
	}else{
    	    system("$ENV{mPUMA}/bin/pick_random_reads.pl",
                '-n',   $number,
                '-r',   $src,
                '-o',   $dest) == 0 or die "Error sub-sampling $src to $number reads and writting to $dest: $!";
	}
}

#############################################################################
# find_minimum - Takes a list of numerical read values and calculates the   #
#               min value.                                                  #
#############################################################################

sub find_minimum {
	my @sortedReads = sort { $a <=> $b } @_;
	return shift @sortedReads;
}

#############################################################################
# find_median - Takes a list of numerical read values and calculates the    #
#               median value.                                               #
#############################################################################
# Input Parameters:                                                         #
#                   @medianReads - The number of reads list to parse        #
# Output Parameters:                                                        #
#                   $median - The median number of reads                    #
#                                                                           #
#############################################################################

sub find_median {
	my @medianReads = @_;
	# Calculate median number of reads.
	@medianReads = sort { $a <=> $b } @medianReads;
	my $median;
	if( (@medianReads % 2) == 1 ) {
		$median = $medianReads[((@medianReads+1)/2)-1];
	} 
	else {
		$median = ($medianReads[(@medianReads/2)-1] + $medianReads[@medianReads/2])/2;
	}
	return $median;
}
	
#############################################################################
# generateDiversityCDHit - Generates temporary diversity files for each MID #
# directory in the projects directory for input to the                      #
# generateDiversityOnlyIsotigs subroutine. Scaled to the median number of   #
# reads in the library.	                                                    #
#############################################################################
# Input Parameters:                                                         #
#                   $assemblyDir - The assembly directory.                  #
#                   $assemblyWateredBLAST - Non-redundant OTU waterBlast.   #
#		    $totalReadSet - "total.fna" format output file.	    #
#                   $totalReadSetWB - The waterBlast "total.fna" format     #
#                                     output file.                          #
#                   $lDir - The library directory.                          #
#                   $scale - The number of reads list to parse.             #
#                   $clusterReport - The cluster-cd-hit report file.        #
#                   $seqcleanFile - Non-redundant set of OTU sequences.     #
# Output Parameters:                                                        #
#                   $output - The DiversityCDHit format output file.        #
#                                                                           #
#############################################################################

sub generateDiversityCDHit{

	my $assemblyDir = shift;
	die "Error lost assemblyDir" unless defined $assemblyDir;
	
	my $assemblyWateredBLAST = shift;
	die "Error lost watered BLAST of the non-redundant OTU" unless defined $assemblyWateredBLAST;
	
	my $totalReadSet = shift;
	die "Error lost total reads" unless defined $totalReadSet;

	my $totalReadSetWB = shift;
	die "Error lost total read wateredBLAST" unless defined $totalReadSetWB;

	my $lDir = shift;
	die "Error lost the directory for this library" unless defined $lDir;
	my $readsIdsOutputFile = join('/',$lDir,$read_ids_file);

	my $scale = shift;
	die "Error lost value to scale to" unless defined $scale;

	my $clusterReport = shift;
	die "Error lost cluster report" unless defined $clusterReport;

	my $seqclean = shift;
	die "Error lost seqclean file" unless defined $seqclean;

	warn "$lDir....\n" if($VERBOSE);

	my $output = join('/',$lDir,'diversity.seqclean-CDhit.scaled.'.$scale.'.txt');

	unless(-s $output){
		system($generateDiversity,
			'-a',		$assemblyDir,
			'-aw',		$assemblyWateredBLAST,
			'-i',		$totalReadSet,
			'-iw',		$totalReadSetWB,
			'-s',		$readsIdsOutputFile,
			'-scale',	$scale,
			'-c',		$clusterReport,
			'-sc',		$seqclean,
			'-v',		0,
			'-o', 		$output
		) == 0 or die "Error calling $generateDiversity for $lDir: $!";
	}
	return $output;
}

sub generateDiversityToc{

	my $lDir = shift;
	die "Error lost the directory for this library" unless defined $lDir;
	
	my $totalToc = shift;
	die "Error lost total TOC" unless defined $totalToc;

	my $assemblyWateredBLAST = shift;
	die "Error lost watered BLAST of the non-redundant OTU" unless defined $assemblyWateredBLAST;
	
	my $clusterReport = shift;

	# What are we reading in order to generate what? 
	my $inputReads = join('/',$lDir,$read_sub_file);
	my $output = join('/',$lDir,$diversity_file);

	# setup the arguments 
	my @args = (
		'-t',	$totalToc,
		'-i',	$inputReads,
		'-w',	$assemblyWateredBLAST,
		'-o',	$output
	);
	# add the clustering report if its defined
	push(@args,'-c',$clusterReport) if defined $clusterReport;
	
	# calculate the diversity if its not already there
	unless(-s $output){
		warn "Generating diversity from $inputReads => $output" if $VERBOSE;
		system($generateDiversity,@args) == 0 or die "Error calling $generateDiversity for $lDir: $!";
	}
	return $output;
}




#############################################################################
# getNUM_CPU - Gets the CPU count on the server and multiplies the total    #
# number of CPUs by a predesignated percentage. (rounds down if CPU count   #
# is a fraction)                                                            #
#############################################################################
# Input Parameters:                                                         #
# 		    $assembly_dir - The assembly directory		    #
# 		    $percent_CPUs - The percent of total CPUs to allocate.  #
# Output Parameters:                                                        #
#                   $NUM_CPU - Percentage of CPUs on the server.            #
#                                                                           #
#############################################################################

sub getNUM_CPU{
	my ($assembly_dir, $getNUM_CPUCmd, $pid, $NUM_CPU, $percent_CPUs);
	$assembly_dir = shift;
        die "Error lost the assembly directory" unless defined $assembly_dir;
	$percent_CPUs = shift;
        die "Error lost the CPU allocation percent" unless defined $percent_CPUs;

	$getNUM_CPUCmd = "cat /proc/cpuinfo";
	if($percent_CPUs =~ m/(\d+.\d+)/){
		$percent_CPUs = $1;
	}
	elsif($percent_CPUs =~ m/\d+/){
		$percent_CPUs = $percent_CPUs/100;
	}
	# getCPUNUM I/O command
	local (*NUM_CPU_OUT, *NUM_CPU_IN);
	$pid = open2(\*NUM_CPU_OUT,\*NUM_CPU_IN, $getNUM_CPUCmd) or die "Error calling open2: $!";
	close NUM_CPU_IN or die "Error closing STDIN to getNUM_CPU process: $!";	
	$NUM_CPU = 2;
	while (<NUM_CPU_OUT>){ 
		if($_ =~ m/^cpu cores\t: (\d+)\n/){
			$NUM_CPU = $1;
			$NUM_CPU = floor(($NUM_CPU + 1)*$percent_CPUs);
			last;
		}
	}
	
	close NUM_CPU_OUT or die "Error closing STDOUT from getNUM_CPU process: $!";
	wait;	
	return $NUM_CPU;
}


#############################################################################
# generateDiversityOnlyIsotigs - Creates the actual diversity files for     #
# each MID directory in the projects directory for downstream analysis.	    #
#############################################################################
# Input Parameters:                                                         #
#                   $lib - The library name. The path to diversity file.    #
#                   $file - The path to diversity file.                     #
# Output Parameters:                                                        #
#                   $output - The DiversityOnlyIsotigs format output file.  #
#                                                                           #
#############################################################################

sub generateDiversityOnlyIsotigs{
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $lib = shift;
	die "Error lost the name of this library" unless defined $lib;

	my ($filename,$base_dirs) = fileparse($file);
	chop($base_dirs);
	my $output = join('/',$base_dirs,$lib . '.txt'); # output is in the same dir but is just named lib.txt. This needs to then get converted for genespring etc but will have the correct name
	warn "$output existed and is being overwritten" if ((-s $output) and ($VERBOSE));

	# read in all the lines starting with isotig
	open(INPUT,"<$file") or die "Error opening $file: $!";
	my @lines;
	while(<INPUT>){
		if ($_ =~ /^isotig/){
			push(@lines,$_);
		}
	}
	close INPUT or die "Error closing $file: $!";

	# write them out in a sorted fashion
	open(OUTPUT,">$output") or die "Error opening $output: $!";
	foreach my $line (sort @lines){
		print OUTPUT $line;
	}
	close OUTPUT or die "Error clsoing $output: $!";

	return $output;
}

#############################################################################
# generateIsotigList - Generate a list of isotigs that exist in each        #
# diversity file                                                            #
#############################################################################
# Input Parameters:                                                         #
#                   $working_dir - The working directory where the          #
#                                  isotig_list file is located.             #
#                   $file_list - The isotig_list files to parse for         #
#                                isotigs.                                   #
# Output Parameters:                                                        #
#                   $isotig_file - The isotigs list output file.            #
#                                                                           #
#############################################################################

sub generateIsotigList{
	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $file_list = shift;
	die "Error lost list of files" unless defined $file_list;

	warn "\nCreating the isotig list for downstream analysis....\n" if($VERBOSE);

	my %isotigs;
	foreach my $file (@{$file_list}){
		open(FILE,"<$file") or die "Error opening $file: $!";
		while(<FILE>){
			chomp $_;
			if ($_ =~ /^(\w+)/){	# this could probably be changed to grab the first \w+ chars - more generic
				$isotigs{$1}++;
			}	
		}
		close FILE or die "Error closing $file: $!";
	}

	my $isotig_file = join('/',$working_dir,'OTU-list');

	warn "WARNING $isotig_file existed and is being overwritten" 
		if ((-s $isotig_file) and ($VERBOSE));

	open(FILE,">$isotig_file") or die "Error opening $isotig_file for writting: $!";
	foreach my $isotig (sort keys %isotigs){
		print FILE $isotig,"\n";
	}
	close FILE or die "Error closing $isotig_file: $!";
	return $isotig_file;
}

#############################################################################
# generateTechFile - Creates template technology file for OTU extraction    #
# used by the convertTech sub routine.                                      #
#############################################################################
# Input Parameters:                                                         #
#                   $isotig_list_file - The isotigs list format file.       #
#                   $scaleValue - The number of median reads.               #
# Output Parameters:                                                        #
#                   $tech_file - The path to the technology file            #
#                                                                           #
#############################################################################

sub generateTechFile{

	my $isotig_list_file = shift;
	die "Error lost isotig file" unless defined $isotig_list_file;
	my $scaleValue = shift;
	die "Error lost the number to scale to" unless defined $scaleValue;

	my $tech_file = $isotig_list_file . '.tech.tsv';
	
	warn "\nCreating the technology.tsv file for OTU extraction....\n" if($VERBOSE);
	unless(-s $tech_file){
		open(INPUT,"<$isotig_list_file") or die "Error opening $isotig_list_file: $!";
		open(OUTPUT,">$tech_file") or die "Error opening $tech_file for writting: $!";
	
		print OUTPUT "OTU\tPercent\tScaled$scaleValue\tActualCount\tLabel\n";	
		while(<INPUT>){
			chomp $_;
			if ($_ =~ /^(\w+)/){
				print OUTPUT "$1\t0\t0\t0\tsomething\n";
			}else{
				die "Error problem with line from $isotig_list_file";
			}
		}
		close INPUT or die "Error closing $isotig_list_file: $!";
		close OUTPUT or die "Error closing $tech_file: $!";
	}
	return $tech_file;	
}

#############################################################################
# convertTech - Converts the template files into technology (*.tsv) files   #
# usable by gene_spring.                                                    #
#############################################################################
# Input Parameters:                                                         #
#                   $working_dir - The working directory where the          #
#                                  output files are going.                  #
#                   $tech_file - The path to the technology file            #
#                   $files2Convert - hashref which maps to the files to     #
#                                    convert                                #
# Output Parameters:                                                        #
#                   $status - Status of the convertTech check               #
#                   $genespring_dir - The genespring directory              #
#                                                                           #
#############################################################################

sub convertTech{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $tech = shift;
	die "Error lost technology file" unless defined $tech;
	
	my $libDir = shift;
	die "Error lost library dir" unless defined $libDir;

	my $libraries = shift;
	die "Error lost libraries" unless defined $libraries;


	# setup the output directory
	my $genespring_dir = join('/',$working_dir,'FILES_FOR_GENESPRING');
	unless(-d $genespring_dir){
		mkdir($genespring_dir, 0777) or die "Can't make directory: $!";
	}

	# read the tech file for OTU order
	my @OTU_order;
	open(TECH,"<$tech") or die "Error opening $tech: $!";
	my $header = <TECH>;
	chomp $header;
	while(<TECH>){
	        chomp $_;
	        my ($otu) = split(/\t/,$_);
	        push(@OTU_order,$otu);
	}
	close TECH or die "Error closing $tech: $!";

     
	foreach my $libPath (@{$libraries}){
                my $lib = basename($libPath);
                my $file = join('/',$libDir,$libPath,$diversity_file);

		# Read in the existing OTU infoarmation for this library
	        open(FILE,"<$file") or die "Error opening $file for reading: $!";
	        my %OTU_data;
	        while(<FILE>){
	                chomp $_;
	                my ($otu) = split(/\t/,$_);
	                $OTU_data{$otu} = $_;
	        }
	        close FILE or die "Error closing $file: $!";
	
	        # write out the OTU information in EXACTLY the same order as the technology file
	        my $outfile = join('/',$genespring_dir,$lib . '.txt');
	        open(DEST,">$outfile") or die "Error opening $outfile for writting: $!";
	        print DEST $header,"\n";
	        foreach my $otu (@OTU_order){
	#               my $line = join("\t",$otu,'N/A','N/A','N/A','something');
	                my $line = join("\t",$otu,0,0,0,'something');
	                if (defined $OTU_data{$otu}){
	                        $line = $OTU_data{$otu};
	                }
	                print DEST $line,"\n";
	        }
	        close DEST or die "Error closing $outfile: $!";
	}

	# Checks the convertedTech output formated file for errors.
	my $status = convertTechCheck($genespring_dir, $tech);
	return ($genespring_dir,$status);
}


#############################################################################
# convertTechCheck - Checks the individual "diversityTech" output files     #
# for errors.                                                               #
# If success returns successful message                                     #
# Otherwise kills program and prints the file/line number of the error.     #
#############################################################################
# Input Parameters:                                                         #
#                   $genespring_dir - The genespring directory              #
#                   $tech_file - The path to the technology file            #
#                                                                           #
# Output Parameters:                                                        #
#                   $status - Status of the convertTech check               #
#                                                                           #
#############################################################################

sub convertTechCheck{

	my $genespring_dir = shift;
	die "error lost directory" unless defined $genespring_dir;
	my $tech_file = shift;
	die "Error lost technolgy file" unless defined $tech_file;

	my $status = -1;

	my @isotig_list;
	open(FILE,"<$tech_file") or die "Error opening $tech_file for reading: $!";
	my $header = <FILE>;
	while(<FILE>){
		chomp $_;
		my ($id) = split(/\t/,$_);
		push(@isotig_list,$id);
	}
	close FILE or die "Error closing $tech_file: $!";
	
	opendir(DIR,$genespring_dir) or die "Error opening $genespring_dir: $!";
	while(my $file = readdir(DIR)){
		next if ($file =~ m/^\./);
		open(FILE,"<$genespring_dir/$file") or die "Error opening $genespring_dir/$file: $!";
		my $head = <FILE>;
		# could check if head eq header 
		my $i = 0;
		while (<FILE>){
			chomp $_;
			my @line = split("\t", $_);
			if(scalar(@line) == 5){
				if(($line[1] =~ m/\d+.\d+|\d|N\/A/) and ($line[2] =~ m/\d+.\d+|\d|N\/A/) and ($line[3] =~ m/\d+.\d+|\d|N\/A/) ){
					$status = 1;
				}
			}else{
				$status = 0;
				if ($status == 0){
					my $error = "Error in file: " . join('/',$genespring_dir,$file) . " line# " . ($i + 1) . " status: format failure\n" . $_;
					$status = $error;
					return $status;
				}
				
			}
			# check the ID order
			if ($isotig_list[$i] ne $line[0]){
				$status = "Error in $tech_file at line " . ($i + 1) . " id mismatch $isotig_list[$i] with technology file $line[0]\n" . $_;
				return $status;
			}
			$i++;
		}
		close FILE or die "Error closing $genespring_dir/$file: $!";
	}
	closedir(DIR) or die "Error closing $genespring_dir: $!";
	return $status;
}

#############################################################################
# generateFastaMEGAN - Creates a MEGAN FASTA input file for each library    #
#############################################################################
# Input Parameters:                                                         #
#                   $working_dir - The working directory to construct a     #
#                                  MEGAN directory where the output files   #
#                                  are going.                               #
#                   $divfile - The path to diversity file.                  #
#                   $CDHitClustFnaFile - The path to the clstr.fna file.    #
#############################################################################

sub generateFastaMEGAN{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $divfile = shift;
	die "Error lost list of files" unless defined $divfile;
	my $CDHitClustFnaFile = shift;
	die "Error lost clstr.fna file" unless defined $CDHitClustFnaFile;

	my $megan_files_dir = join('/',$working_dir,'FILES_FOR_MEGAN');
	unless(-d $megan_files_dir){
		mkdir($megan_files_dir, 0777) or die "Can't make directory: $!";
	}

	my @libfile = split('/', $divfile);
	my $libName;
	if ( $libfile[-2] =~ m/(.+)/){	# this could probably be changed to grab the first \w+ chars - more generic
		$libName = $1;
	}		

	# convert the files
	my $meganLibDir = join( '/', $megan_files_dir, $libName);
	my $meganFile = $meganLibDir . '.megan.fasta';
	my $createMeganInputCmd = "$createMeganInput -d $divfile -f $CDHitClustFnaFile -o $meganFile\n\n";
	warn $createMeganInputCmd if $VERBOSE;

	unless(-s $meganFile){
		system($createMeganInput,
			'-d', $divfile,
			'-f', $CDHitClustFnaFile,
			'-o', $meganFile
		) == 0 or die "Error calling $createMeganInput: $?";
	}
}

#############################################################################
# generateBLASTX -  Creates a BLASTx input file for the generateBlastxMEGAN #
# and translateStranded sub-routines.                                       #
#############################################################################
# Input Parameters:                                                         #
#                   $FnaFile - The path to the fna file.                    #
# Output Parameters:                                                        #
#                   $FnaBlastxFile - The path to the                        #
#                                     fna.blastx file.                      #
#############################################################################

sub generateBLASTX{
	my $FnaFile = shift;
	die "Error lost fna file" unless defined $FnaFile;
	my $vOption = shift;
	die "Error lost vOption" unless defined $vOption;
	my $bOption = shift;
	die "Error lost bOption" unless defined $bOption;
	my $FnaBlastxFile = $FnaFile . '.blastx';

	formatdbPep($mPUMA::cpnDB_nr_pep);

	# convert the files
	unless(-s $FnaBlastxFile){	
		warn "\nGenerating BLASTx file....\n" if($VERBOSE);
		my $blastallCmd  = "blastall -p blastx -i $FnaFile -d $mPUMA::cpnDB_nr_pep -F F -v $vOption -b $bOption -o $FnaBlastxFile -a $cpu_seqclean\n\n";
		warn $blastallCmd if $VERBOSE;
		my $status = system($blastall, 
		'-p', 'blastx', 
		'-i', $FnaFile, 
		'-d', $mPUMA::cpnDB_nr_pep, 
		'-F', 'F', 
		'-v', $vOption, 
		'-b', $bOption, 
		'-o', $FnaBlastxFile,
		'-a', $cpu_cores
		);
		if ($status ne 0){
			if (-e $FnaBlastxFile){
				warn "Error calling $blastall: $?";
				warn "$FnaBlastxFile exists but could have partial results so trying to delete it";
				unlink($FnaBlastxFile) or die "Error removing $FnaBlastxFile: $!";
				die "Aborting pipeline because of problem with $blastall: $?";
			}
                }


	}
	return $FnaBlastxFile;
}

#############################################################################
# generateBlastxMEGAN - Creates a MEGAN BLASTX input file for each library  #
#############################################################################
# Input Parameters:                                                         #
#                   $working_dir - The working directory to construct a     #
#                                  MEGAN directory where the output files   #
#                                  are going.                               #
#                   $divfile - The path to diversity file.                  #
#                   $CDHitClustFnaBlastxFile - The path to the              #
#                                              clstr.fna.blastx file.       #
#############################################################################

sub generateBlastxMEGAN{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $divfile = shift;
	die "Error lost list of files" unless defined $divfile;
	my $CDHitClustFnaBlastxFile = shift;
	die "Error lost clstr.fna.blastx file" unless defined $CDHitClustFnaBlastxFile;

	my $megan_files_dir = join('/',$working_dir,'FILES_FOR_MEGAN');
	unless(-d $megan_files_dir){
		mkdir($megan_files_dir, 0777) or die "Can't make directory: $!";
	}

	my @libfile = split('/', $divfile);
	my $libName;
	if ( $libfile[-2] =~ m/(.+)/){	# this could probably be changed to grab the first \w+ chars - more generic
		$libName = $1;	
	}		

	# convert the files
	my $meganLibDir = join( '/', $megan_files_dir, $libName);
	my $meganFile = $meganLibDir . '.megan.blastx';
	my $createMeganBlastxInputCmd = "$createMeganInput -d $divfile -b $CDHitClustFnaBlastxFile -o $meganFile\n\n";
	warn $createMeganBlastxInputCmd if $VERBOSE;

	unless(-s $meganFile){
		system($createMeganBlastxInput,
			'-d', $divfile,
			'-b', $CDHitClustFnaBlastxFile,
			'-o', $meganFile
		) == 0 or die "Error calling $createMeganBlastxInput: $?";
	}
}

sub generateCSVInputMEGAN{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $libName = shift;
	die "Error lost library name " unless defined $libName;
	my $divfile = shift;
	die "Error lost list of files" unless defined $divfile;
	my $CDHitClustFnaBlastxFile = shift;
	die "Error lost clstr.fna.blastx file" unless defined $CDHitClustFnaBlastxFile;

	my $megan_files_dir = join('/',$working_dir,'FILES_FOR_MEGAN');
	unless(-d $megan_files_dir){
		mkdir($megan_files_dir, 0777) or die "Can't make directory: $!";
	}

	my $possible_megan_results = join( '/', $megan_files_dir, $libName . '.megan');
	if (-s $possible_megan_results){
		warn "Warning it appears that $possible_megan_results exists so skipping";
		return $possible_megan_results;
	}

	# convert the files
	my $meganFile = join( '/', $megan_files_dir, $libName . '.gz');

	my $createMeganInputCmd = "$createMeganCSVInput -d $divfile -b $CDHitClustFnaBlastxFile -o $meganFile\n\n";
	warn $createMeganInputCmd if $VERBOSE;
	
	unless(-s $meganFile){
		system($createMeganCSVInput,
			'-d', $divfile,
			'-b', $CDHitClustFnaBlastxFile,
			'-o', $meganFile
		) == 0 or die "Error calling $createMeganCSVInput: $?";
	}
	return $meganFile;
}


sub generateClassifierInputMEGAN{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $libName = shift;
	die "Error lost library name " unless defined $libName;
	my $divfile = shift;
	die "Error lost list of files" unless defined $divfile;
	my $classifierFile = shift;
	die "Error lost classifier file" unless defined $classifierFile;

	# ensure that the MEGAN directory exists
	my $megan_files_dir = join('/',$working_dir,'FILES_FOR_MEGAN');
	unless(-d $megan_files_dir){
		mkdir($megan_files_dir, 0777) or die "Can't make directory: $!";
	}

	my $possible_megan_results = join( '/', $megan_files_dir, $libName . '.megan');
	if (-s $possible_megan_results){
		warn "Warning it appears that $possible_megan_results exists so skipping";
		return $possible_megan_results;
	}

	# convert the files
	my $meganFile = join( '/', $megan_files_dir, $libName);

	unless(-s $meganFile){
		warn "Generating MEGAN input using the classifier results and diversity info from $divfile" if $VERBOSE;
		system($createMeganClassifierInput,
			'-d', $divfile,
			'-f', $classifierFile,
			'-o', $meganFile
		) == 0 or die "Error calling $createMeganClassifierInput: $?";
	}
	return $meganFile;
}

sub generateClassifierProfilesFromMegan{

        my $working_dir		= shift;
        die "Error lost directory" unless defined $working_dir;
        my $libName		= shift;
        die "Error lost library name " unless defined $libName;
        my $meganClassifierFile	= shift;
        die "Error lost classifier file" unless defined $meganClassifierFile;
	my $taxonomicTerm	= shift;
	$taxonomicTerm		= 'phylum' unless defined $taxonomicTerm;
	my $percentConfidence	= 0.80;	# default is defined here because this is the place where it should be set based on some kind of established criteria	

        # ensure that the MEGAN directory exists
        my $output_dir = join('/',$working_dir,'CLASSIFIER_PROFILES');
        unless(-d $output_dir){
                mkdir($output_dir, 0777) or die "Can't make directory: $!";
        }

	my $output_file = join('/',$output_dir,$libName .'.'.$taxonomicTerm);
	unless (-s $output_file){
		system("$ENV{mPUMA}/bin/convert_classifier_to_profile.pl",
			'-c',	$meganClassifierFile,
			'-t',	$taxonomicTerm,
			'-p',	$percentConfidence,
			'-o',	$output_file) == 0 or die "Error calling convert_classifier_to_profile.pl -c $meganClassifierFile -t $taxonomicTerm -p $percentConfidence -o $output_file: $!";
	}
	return $output_file;
}

# This is a routine to analyze the CSV Megan input files
sub convertMEGANCSV{
	my $file = shift;
	die "Error lost file to try and convert" unless defined $file;
	my $unlink = shift;
	$unlink = 0 unless defined $unlink;

	if (defined $ENV{DISPLAY}){
		# there should be support for X and so this should work 
		my $output = $file . '.megan';
		if (-s $output){
			warn "Warning skipping MEGAN analysis of $file because $output already exists";
		}else{
			my ($filename,$dir) = fileparse($file);
			my $outfilename = $filename . '.megan';
			my $cmd = "MEGAN +g false +s false -S true +w false -x \"set dir=$dir;import csv=reads separator=comma file=$filename taxonomy=true toppercent=10.0 minscore=0 minsupport=5; save file=$outfilename summary=true; quit;\"";
			warn $cmd if $VERBOSE;
			system($cmd) == 0 or die "Error calling $cmd: $!";
			if ((-s $output)and($unlink)){
				unlink($file) or die "Error removing $file: $!";
			}
				# output was created so we can remove the input file
		}
		return $output;
	}else{
		warn "NON-fatal error: you are trying to use megan in a command-line fashion and it needs a terminal which has X11 support. You will need to call this from a terminal which has a valid DISPLAY environment variable set: skipping $file";
	}
	return $file;
}

#############################################################################
# generateInputMOTHUR - Creates a MOTHUR input file consisting of each      #
# library                                                                   #
#############################################################################
# Input Parameters:                                                         #
#                   $working_dir - The working directory to construct a     #
#                                  MOTHUR directory where the output file   #
#                                  is going.                                #
#                   $genespring_dir - The genespring directory              #
#############################################################################

sub generateInputMOTHUR{
	
	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	
	my $libDir = shift;
	die "Error lost library dir" unless defined $libDir;

	my $libraries = shift;
	die "Error lost libraries" unless defined $libraries;

	my $mothur_files_dir = join('/',$working_dir,'FILES_FOR_MOTHUR');
	unless(-d $mothur_files_dir){
		mkdir($mothur_files_dir, 0777) or die "Can't make directory: $!";
	}	
	my $mothur_file = join('/',$mothur_files_dir,'freq_mothur.txt');
	
	my %OTUabundance;
	my %OTUnames;
	
	foreach my $libPath (@{$libraries}){
		my $lib = basename($libPath);
		my $divFile = join('/',$libDir,$libPath,$diversity_file);
		print "\t$lib => $divFile\n";
		     
		open(FILE,"<$divFile") or die "Error opening $divFile for reading: $!";
        
		while(<FILE>){
			chomp $_;
			my @parts = split(/\s+/,$_);
			my $OTU = $parts[0];
			my $actualCount = $parts[3];
			next unless $actualCount >= 1;
			$OTUabundance{$lib}->{$OTU} = $actualCount;
			$OTUnames{$OTU}++;
		}
		close FILE or die "Error closing $divFile: $!";
	}

	# write out the output
	open(OUTPUT,">$mothur_file") or die "Error opening $mothur_file for writting: $!";
	my @sortedOTUnames = sort {return $a cmp $b} keys %OTUnames; # sort the OTUnames so that we can ensure that the columns in the MOTHUR input file will line up accordingly
	my $numberOfOTU = @sortedOTUnames;

	# I wish mothur would allow some type of commenting in the file format so there there could be some reuse of this datafile...
	#print OUTPUT "#",join("\t",'NULL','LIBRARY NAME','#OTUs',@sortedOTUnames),"\n";
	foreach my $lib (keys %OTUabundance){
	
                # create the ordered set of abundance values
                my @orderedAbundances;
                foreach my $id (@sortedOTUnames){
                        my $abundance = 0;      # if there was no abundance (or it was below $min) then set this to 0
                        $abundance = $OTUabundance{$lib}->{$id} if defined $OTUabundance{$lib}->{$id};
			warn "$lib,$id,$abundance";
                        push(@orderedAbundances,$abundance);
                }

                # ensure that there are the same number of abundances as we expect!!!
                my $numAbundances = @orderedAbundances;
                die "Error $numAbundances != $numberOfOTU and so this would corrupt MOTHUR input file"
                        unless ($numAbundances == $numberOfOTU);

                print OUTPUT join("\t",'NA',$lib,$numberOfOTU,@orderedAbundances),"\n";
        }
	
	close OUTPUT or die "Error closing $mothur_file: $!";

	return ($mothur_file,$mothur_files_dir);

}
sub generateInputMOTHURold{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $genespring_dir = shift;
	die "Error lost genespring directory" unless defined $genespring_dir;

	my $mothur_files_dir = join('/',$working_dir,'FILES_FOR_MOTHUR');
	unless(-d $mothur_files_dir){
		mkdir($mothur_files_dir, 0777) or die "Can't make directory: $!";
	}	
	my $mothur_file = join('/',$mothur_files_dir,'freq_mothur.txt');
	# convert the files
	my $createMothurInputCmd = "$createMothurInput -d $genespring_dir -o $mothur_file\n\n";
	warn $createMothurInputCmd if $VERBOSE;

	unless(-s $mothur_file){
		system($createMothurInput,
			'-d', $genespring_dir,
			'-o', $mothur_file,
			'-m', 2, # this is the minimum # of times an OTU must be seen in order to be included in the MOTHUR input file(s)
		) == 0 or die "Error calling $createMothurInput: $?";
	}
	return ($mothur_file, $mothur_files_dir);
}

#############################################################################
# processInputMOTHUR - Processes MOTHUR input files that consist of each    #
# library performs calc=nseqs-coverage-npshannon-simpson-sobs-chao and      #
# rarefaction data                                                          #
#############################################################################
# Input Parameters:                                                         #
#                   $mothur_file - The working directory to construct a     #
#                                  MOTHUR directory where the output file   #
#                                  is going.                                #
#                   $mothur_file_dir - The genespring directory             #
#############################################################################

sub processInputMOTHUR{

	my $mothur_file = shift;
	die "Error lost MOTHUR input file" unless defined $mothur_file;
	my $mothur_files_dir = shift;
	die "Error lost MOTHUR files directory" unless defined $mothur_files_dir;

	my $mothur_outfiles_dir = join('/',$mothur_files_dir,'MOTHUR_OUTPUT_FILES');
	if (-d $mothur_outfiles_dir){
		warn "$mothur_outfiles_dir exists so assuming that mothur has been run already" if $VERBOSE;
	}else{
		mkdir($mothur_outfiles_dir, 0777) or die "Can't make directory: $!";
		my $processInputMothurCmd = "$mothur \"#set.dir(output=$mothur_outfiles_dir);summary.single(shared=$mothur_file, calc=nseqs-coverage-npshannon-simpson-sobs-chao);rarefaction.single(shared=$mothur_file)\"";
		warn "$processInputMothurCmd\n\n" if $VERBOSE;
		system($processInputMothurCmd) == 0 or die "Error calling $processInputMothurCmd: $?";
	}
}

#############################################################################
# generateInputUNIFRAC - Creates a UNIFRAC input file consisting of each    #
# library                                                                   #
#############################################################################
# Input Parameters:                                                         #
#                   $working_dir - The working directory to construct a     #
#                                  UNIFRAC directory where the output file  #
#                                  is going.                                #
#                   $genespring_dir - The genespring directory              #
#############################################################################
sub generateInputUNIFRAC{

	my $working_dir = shift;
	die "Error lost directory" unless defined $working_dir;
	my $diversityFiles = shift;
	die "Error lost diversityFiles" unless defined $diversityFiles;
	my $fastaFile = shift; 
	die "Error lost fasta" unless defined $fastaFile;
	my $seqType = shift;
	die "Error lost sequence type" unless defined $seqType;
	die "Error sequence type $seqType is not supported" unless (($seqType eq 'Nucleotide')or($seqType eq 'Protein'));

	my $unifrac_files_dir = join('/',$working_dir,'FILES_FOR_UNIFRAC');
	unless(-d $unifrac_files_dir){
		mkdir($unifrac_files_dir, 0777) or die "Can't make directory: $!";
	}

	# generate a quick multiple sequence alignment
	my $multipleSeqAln = t_coffee_fasta($fastaFile);

	# do a phylogenetic workup using FastTree
	my $fastTree;
	if ($seqType eq 'Nucleotide'){
		$fastTree = fastTree_nuc($multipleSeqAln,1); # call FastTree using gtr
	}elsif($seqType eq 'Protein'){
		$fastTree = fastTree_aln($multipleSeqAln); # call FastTree
	}else{
		die "Unexpected case where seqtype = $seqType";
	}

	# symlink the phylogenetic tree into the Unifrac dir
	my $unifracTree = join('/',$unifrac_files_dir,'FastTree');
	unless (-e $unifracTree){
		symlink($fastTree,$unifracTree) or die "Error creating symlink from $fastTree => $unifracTree: $!";
	}

	# create a hash which goes from library Name to diversity file
	my %library2diversityMap;
	foreach my $f (@{$diversityFiles}){
		my @parts = split(/\//,$f);
		my $file = pop @parts;
		if ($file eq $diversity_file){
			my $lib = pop @parts;
			die "Error multiple diversity files defined for $lib" if (defined $library2diversityMap{$lib});
			$library2diversityMap{$lib} = $f;
		}else{
			die "Error last part of $f was [$file]. Was expecting $diversity_file";
		}
	}		

	# setup a simple Category file for UniFrac 
	my $unifracCategoryFile = join('/',$unifrac_files_dir,'Categories.txt');
	unless (-s $unifracCategoryFile){
		open(FILE,">$unifracCategoryFile") or die "Error opening $unifracCategoryFile for writting: $!";
		print FILE "#SampleID	Subcategory1	Description\n";
		print FILE "# The user MUST enter categorical information as per the UniFrac documentation http://bmf2.colorado.edu/fastunifrac/tutorial.psp\n";
		print FILE "# The information here is just dummy placeholder info.\n";
		foreach my $lib (sort {return $a cmp $b} keys %library2diversityMap){
			print FILE join("\t",$lib,'SubCategory1',"Description of $lib goes here"),"\n";
		}
		close FILE or die "Error closing $unifracCategoryFile: $!";
	}

	# create the sample ID mapping file from the diversity files
	my $unifracSampleIDfile = join('/',$unifrac_files_dir,'sampleIDmappingFile.txt');
	unless(-s $unifracSampleIDfile){
		my %sampleToLibCount;
		foreach my $lib (keys %library2diversityMap){
			open(DIV,"<$library2diversityMap{$lib}") or die "Error opening $library2diversityMap{$lib} for reading: $!";
			while(<DIV>){
				chomp $_;
				my ($OTU,$percent,$scaled,$actual,$match) = split(/\t/,$_);
				die "Error got multiple counts for $OTU from $lib" if (defined $sampleToLibCount{$OTU}->{$lib});
				$sampleToLibCount{$OTU}->{$lib} = $actual;
			}
			close DIV or die "Error closing $library2diversityMap{$lib}: $!";
		}
		open(SAMPLE,">$unifracSampleIDfile") or die "Error opening $unifracSampleIDfile for writting: $!";
		foreach my $OTU (sort {return $a cmp $b} keys %sampleToLibCount){
			# this is not going to write a 0 for those libraries which did not observe a specific OTU. This is as per the UniFrac documentation http://bmf2.colorado.edu/fastunifrac/tutorial.psp
			foreach my $lib (sort {return $a cmp $b} keys %{$sampleToLibCount{$OTU}}){
				print SAMPLE join("\t",$OTU,$lib,$sampleToLibCount{$OTU}->{$lib}),"\n";
			}
		}
		close SAMPLE or die "Error closing $unifracSampleIDfile: $!";
	}
	return $unifrac_files_dir;
	die "Uni dir = $unifrac_files_dir. Fasta $fastaFile. Category file $unifracCategoryFile ";
}

sub prepFastqBasedAssembly {
	my $topLevelDir = shift;
        die "Error lost the Trinity directory" unless defined $topLevelDir;
	if (!(-e $topLevelDir)){
		mkdir($topLevelDir, 0777) or die "Can't make directory $topLevelDir: $!";
	}

	# Creates the assembly directory
	# N.B. this is attempting to emulate the gsAssembler structure to keep things reasonably similar across methods
	my $assemblyDir = join('/',$topLevelDir,'assembly');
	if (!(-e $assemblyDir)){
		mkdir($assemblyDir, 0777) or die "Can't make directory $assemblyDir: $!";
	}

	# Adds the sequence runs from the SFF files.
	my $NGSdir = get_compressed_fastq_folder($topLevelDir);

	if (!(-d $NGSdir)){
		mkdir($NGSdir, 0777) or die "Can't make directory $NGSdir: $!";
	}
	my %seenNGS;
#	foreach my $file (@_){
#		my $basename = basename($file);
#		die "Error multiple files with the same basename ($basename)" if (defined $seenNGS{$basename});
#		$seenNGS{$basename}++;
#		my $target = join('/',$NGSdir,$basename.'.gz');
#		if (-e $target){
#			warn "$target already exists so moving on" if $VERBOSE;
#		}else{
#			# need to get the data for $file into $target as FASTQ gzip compressed
#			# N.B these files are compressed but the total output in the assembly directory will not be
#			if ($basename =~ /.*\.sff$/){
#				# need to convert SFF to FASTQ compressed
#				system("$ENV{mPUMA}/bin/sffToFastq.pl",'-i',$file,'-o',join('/',$NGSdir,$basename),'-z',1) == 0 
#					or die "Error converting sff file to fastq: $file: $!";
#			}elsif(($basename =~ /.*\.fq$/i)or($basename =~ /.*\.fastq$/i)){
#				# file is FASTQ but not compressed
#				system('gzip',$file) == 0 or die "Error compressing $file: $!";
#				my $newFile = $file . '.gz';
#				# create link to compressed version of the file
#				symlink($newFile,$target) or die "Error creating link to $basename. Do you have multiple files which would have the same basename? That's not allowed: $!";
#			}elsif(($basename =~ /.*\.fq\.gz$/i)or($basename =~ /.*\.fastq\.gz$/i)){
#				# Already in fastq format so just create a link
#				symlink($file,$target) or die "Error creating link to $basename. Do you have multiple files which would have the same basename? That's not allowed: $!";
#			}else{
#				die "Error do not know how to handle $file";
#			}
#		}
#	}

	foreach my $file (@_){
		my $basename = basename($file);
#		die "Error multiple files with the same basename ($basename)" if (defined $seenNGS{$basename});
		warn "Error multiple files with the same basename ($basename)" if (defined $seenNGS{$basename});
		$seenNGS{$basename}++;

		# need to get the data for $file into $target as FASTQ gzip compressed
		# N.B these files are compressed but the total output in the assembly directory will not be
		if ($basename =~ /.*\.sff$/){
			# need to convert SFF to FASTQ compressed
			system("$ENV{mPUMA}/bin/sffToFastq.pl",'-i',$file,'-o',join('/',$NGSdir,$basename),'-z',1) == 0 
				or die "Error converting sff file to fastq: $file: $!";
		}elsif(($basename =~ /.*\.fq$/i)or($basename =~ /.*\.fastq$/i)){
			# file is FASTQ but not compressed
			system('gzip',$file) == 0 or die "Error compressing $file: $!";
			my $newFile = $file . '.gz';
			# create link to compressed version of the file

			my $target = join('/',$NGSdir,$basename.'.gz');
			if (-e $target){
				warn "$target already exists so moving on" if $VERBOSE;
				next;
			}

			symlink($newFile,$target) or die "Error creating link to $basename. Do you have multiple files which would have the same basename? That's not allowed: $!";
		}elsif(($basename =~ /.*\.fq\.gz$/i)or($basename =~ /.*\.fastq\.gz$/i)){

			my $target = join('/',$NGSdir,$basename);
			if (-e $target){
				warn "$target already exists so moving on" if $VERBOSE;
				next;
			}

			# Already in fastq format so just create a link
			symlink($file,$target) or die "Error creating link to $basename. Do you have multiple files which would have the same basename? That's not allowed: $!";
		}else{
			die "Error do not know how to handle $file";
		}


	}

	return ($assemblyDir,$NGSdir,\%seenNGS);
}





sub runTrinity {
	
	my $trinityDir = shift;
        die "Error lost the Trinity directory" unless defined $trinityDir;

	# make sure all the FASTQ files are in place
	my ($assemblyDir,$NGSdir,$seenNGS) = prepFastqBasedAssembly($trinityDir,@_);

	# make the total input for Trinity...
	my $totalInput = join('/',$assemblyDir,'total.fq');
	if (-e $totalInput){
		warn "$totalInput already exists so moving on" if $VERBOSE;
	}else{
		open(TOTAL,">$totalInput") or die "Error opening $totalInput for writting: $!";
		foreach my $basename(keys %{$seenNGS}){
			my $input;
			if ($basename =~ /\.gz$/){
				$input = join('/',$NGSdir,$basename );
			}else{	
				$input = join('/',$NGSdir,$basename . '.gz');
			}
			my $in = IO::Uncompress::Gunzip->new($input) or die "Error opening $input for reading: $!";
			print TOTAL while(<$in>);
			close $in or die "Error closing $input: $!";
		}
		close TOTAL or die "Error closing $totalInput: $!";
	}

	# try and pre-empt Trinity creating the total input file. This is done because it appears that Trinity may actually be trying to add a /W notation to reads which are 454 reads converted to fastq and then submitted through Trinity... this can cause all sorts of issues as the IDs in the sam files are not going to line up with these. 
	# So the proposed solution is to create the single.fa FASTA input for trinity by hand and then link this to total.fna so the rest of the pipeline will continue to function. 
	# this needs to be retested with a full assembly from uncompressed 454 conterted to Fastq

	if (!(-e join('/',$assemblyDir,'single.fa'))){
		# input for Trinity doesn't exist so create it.
		convertFastqToFasta(join('/',$assemblyDir,'total.fq'),join('/',$assemblyDir,'single.fa'));
	}
	# make a link to the total FASTA input to trinity as this is needed for the post-assembly cleanup and its in single.fa
	if (!(-e join('/',$assemblyDir,'total.fna'))){
#		link(join('/',$assemblyDir,'single.fa'),join('/',$assemblyDir,'total.fna')) or die "Error creating a link to single.fa: $!";
		symlink(join('/',$assemblyDir,'single.fa'),join('/',$assemblyDir,'total.fna')) or die "Error creating a link to single.fa: $!";
	}

	# run Trinity
	my $trinityOutput = join('/',$assemblyDir,'Trinity.fasta');
	if (!(-e $trinityOutput)){
		system($trinity,
			'--seqType',	'fq',
			'--JM',		'20G',
			'--single',	$totalInput,
			'--CPU',	$cpu_cores,
			'--output',	$assemblyDir) == 0
			or die "Error calling Trinity: $!";
		# N.B. an untested option would be to attempt to use paired end data by ensuring that the FASTQ files had the proper /1 /2 designations and then telling Trinity to use --run_as_paired here. I think that will work. This will likely be tested shortly using MiSeq data.  
	}else{
		warn "Trinity appears to have already completed so skipping it " if $VERBOSE;
	}	
#	# make a link to the total FASTA input to trinity as this is needed for the post-assembly cleanup and its in single.fa
#	if (-e join('/',$assemblyDir,'single.fa')){
#		if (!(-e join('/',$assemblyDir,'total.fna'))){
#			link(join('/',$assemblyDir,'single.fa'),join('/',$assemblyDir,'total.fna')) or die "Error creating a link to single.fa: $!";
#		}
#	}

	return $assemblyDir;
}

sub convertFastqToFasta {
	my $srcFastq = shift;
	die "Error lost source Fastq" unless defined $srcFastq;
	my $destFasta = shift;
	die "Error lost destination Fasta" unless defined $destFasta;
	my $format = shift;
	$format = 'fastq-illumina' unless defined $format; #allows override for the format of Fastq

	my $inio	= Bio::SeqIO->new(-file => "<$srcFastq", -format => $format);
	my $outio	= Bio::SeqIO->new(-file => ">$destFasta", -format => 'fasta');

	while(my $seq = $inio->next_seq()){
		$outio->write_seq($seq);
	}
}



=head1 mPUMA.pm

Microbial Profiling Using Metagenomic Assembly - Perl Module

=head1 DESCRIPTION

mPUMA.pm - The microbial Profiling Using Metagenomic Assembly perl module consists of a sub-routine tool kit used to construct general or specific purpose pipelines for metagenomic analysis of CPN60 454 data.

=head1 SYNOPSIS

B<mPUMA.pm> - L<Sub-Routine List>

=head2 B<parse_library_file()> - Simple parser which reads in a 2-column file which defines which SFFs are for which libraries and is tab delimited. E.g. LibA\t/path/to/SFF1.sff, LibB\t/path/to/SFF2.sff, LibB\t/path/to/SFF3.sff.

I<Input Parameters>:

B<$libfile>  - The library file for parsing SFFs

I<Output Parameters>:

B<\%info> - Hash reference pointing to the SFF filename paths


=head2 B<newAssembly()> - Creates the assembly directory and populates it with files from the cdna library.

I<Input Parameters>:

B<$newblerDir> - The Newbler project directory

I<Output Parameters>:

B<$assemblyDir> - The assembly directory 


=head2 B<addRun()> - Adds sequence runs from the corresponding SFF files.

I<Input Parameters>:

B<$newblerDir> - The Newbler project directory


=head2 B<runProject()> - Initiates the project assembly using newbler's runProject command.

I<Input Parameters>:

B<$assemblyDir> - The assembly directory 

=head2 B<runProjectCheck()> - Checks the outcome of the project assembly for any errors and troubleshoots any issues associated with the error.

I<Input Parameters>:

B<$assemblyDir> - The assembly directory

=head2 B<runProjectError()> - Reports an error associated with the runProject assembly.

=head2 B<runAssembly()> - Initiates the project assembly using newbler's newAssembly, addRun, and runProject commands.

I<Input Parameters>:

B<$newblerDir> - The Newbler project directory

I<Output Parameters>:

B<$assemblyDir> - The assembly directory

=head2 B<post_assembly_cleanup()> - Runs seqclean and cd_hit clustering on the assembled results.

I<Input Parameters>:

B<$assemblyDir> - The assembly directory

I<Output Parameters>:

B<$seqcleanCDHitClustFnaFile> - The cluster-cd-hit.fna file.

B<$seqcleanCDHitClustReportFile> - The cluster-cd-hit report file.

B<$seqcleanFile> - Non-redundant set of OTU sequences

=head2 B<seqcleanFile()> - Calls seqclean on a Fasta file to remove primers

I<Input Parameters>:

B<$fasta> - The Fasta file to remove primers

I<Output Parameters>:

B<$seqcleanFile> - Non-redundant set of OTU sequences

B<$primers> - The primer MSF file to lookup primers to remove.

=head2 B<cluster_cd_hit()> - Performs the cd_hit clustering on the fasta file to obtain the non-redundant OTU seqeunces file.

I<Input Parameters>:

B<$fasta> - The Fasta file to remove primers

I<Output Parameters>:

B<$seqcleanCDHitClustFnaFile> - The cluster-cd-hit.fna file.

B<$seqcleanCDHitClustReportFile> - The cluster-cd-hit report file.

=head2 B<waterBlastIT()> - Runs a waterBlast on a query and target DNA or RNA ( BLASTN ) sequences to obtain a similarity comparison of the sequencies.

I<Input Parameters>:

B<$wbQuery> - The 454Isotigs.fna query file

B<$wbTarget> - The cpndb_nr database target file

B<blast_prog> - must be either blastp or blastn.

I<Output Parameters>:

B<$wbSymlinkFile> - waterBlast format output symbolic link file

B<$status> - Status of the waterBlast check 


=head2 B<WaterBlastCheck()> - Checks the waterBlast output file for errors If success returns a sucessful message Otherwise kills program and prints the file/line number of the error.

I<Input Parameters>:

B<$wbOutputFileName> - waterBlast format output file to check

I<Output Parameters>:

B<$status> - Status of the waterBlast check


=head2 B<solveOrientation()> - Generates the stranded fasta file and returns the path to that file.

I<Input Parameters>:

B<$fasta> - The Fasta file to remove primers

B<$wbOutput> - waterBlast format output file

I<Output Parameters>:

B<$stranded> - stranded fasta format output file

=head2 B<sffinfoFna()> - Creates a total read set file "total.fna" from the SFF files.

I<Input Parameters>:

B<$dir> - The directory to store the "total.fna" file

I<Output Parameters>:

B<$totalReadSet> - "total.fna" format output file

B<$status> - Status of the sffinfoFna check

=head2 B<sffsToFna()> - Creates a total read set file "total.fna" from the SFF files. Calls sffinfo on all of the SFF files and writes the result to $output

I<Input Parameters>:

B<$output> - Name and filepath of the "total.fna" output file.

I<Output Parameters>:

B<$output> - Name of the "total.fna" output file.

=head2 B<sffinfoFnaCheck()> - Checks the "total.fna" output file for errors. If success returns successful message Otherwise kills program and prints the file/line number of the error.

I<Input Parameters>:

B<$file> -  "total.fna" file to check

I<Output Parameters>:

B<$status> - Status of the sffinfoFna check

=head2 B<calcSpSnMetrics()> - Computes the Specificity and Sensitivity on a query and target sequence to obtain metrics validating the waterBlast hits.

I<Input Parameters>:

B<$totalReadSet> - The total.fna query file

B<$contigFileName> - The 454Isotigs.fna query file

B<$assemblyDir> - The assembly directory

I<Output Parameters>:

B<$snspInputFile> - Input for the calSpSnMetricsCheck sub-routine (unaltered "Specificity_Sensitivity" file)

B<$snspOutputFile> - Input for the calSpSnMetricsCheck sub-routine  (empty "Specificity_Sensitivity" file for extracting all useful information separately.)

B<$status> - Status of the calcSpSnMetrics check

=head2 B<calcSpSnMetricsCheck()> - Checks the "Specificity_Sensitivity" output file for errors. If success returns successful message Otherwise kills program and prints the file/line number of the error.

I<Input Parameters>:

B<$snspInputFile> - Input for the calSpSnMetricsCheck sub-routine. (unaltered "Specificity_Sensitivity" file)

B<$snspOutputFile> - Input for the calSpSnMetricsCheck sub-routine. (empty "Specificity_Sensitivity" file for extracting all useful information separately will exclude all N/A lines)

I<Output Parameters>:

B<$status> - Status of the calcSpSnMetrics check

=head2 B<createLibraries()> - Creates a libraries dir and then a sub folder for each library.

I<Input Parameters>:

B<$dir> - The working directory

B<$lib_info>  - Hash reference pointing to the SFF libraries.

I<Output Parameters>:

B<$libraries_dir> - The libraries directory

=head2 B<generateReadsIds()> - Creates the "reads.ids" files for each subsequent SFF directory in the library directory.

I<Input Parameters>:

B<$libDir> - The library directory

B<$sffFiles> - The SFF filename path

I<Output Parameters>:

B<$numReads> - The number of reads in the "reads.ids" file

=head2 B<find_median()> - Takes a list of numerical read values and calculates the median value.

I<Input Parameters>:

B<@medianReads> - The number of reads list to parse

I<Output Parameters>:

B<$median> - The median number of reads


=head2 B<generateDiversityCDHit()> - Generates temporary diversity files for each MID directory in the projects directory for input to the generateDiversityOnlyIsotigs subroutine. Scaled to the median number of reads in the library.

I<Input Parameters>:

B<$assemblyDir> - The assembly directory.

B<$assemblyWateredBLAST> - Non-redundant OTU waterBlast.

B<$totalReadSet> - "total.fna" format output file.

B<$totalReadSetWB> - The waterBlast "total.fna" format output file.

B<$lDir> - The library directory.

B<$scale> - The number of reads list to parse.

B<$clusterReport> - The cluster-cd-hit report file.

B<$seqcleanFile> - Non-redundant set of OTU sequences.

I<Output Parameters>:

B<$output> - The DiversityCDHit format output file.


=head2 B<getNUM_CPU()> - Gets the CPU count on the server and multiplies the total number of CPUs by a predesignated percentage. (rounds down if CPU count is a fraction)

I<Input Parameters>:

B<$assembly_dir> - The assembly directory

B<$percent_CPUs> - The percent of total CPUs to allocate.

I<Output Parameters>:

B<$NUM_CPU> - Percentage of CPUs on the server.


=head2 B<generateDiversityOnlyIsotigs()> - Creates the actual diversity files for each MID directory in the projects directory for downstream analysis.

I<Input Parameters>:

B<$lib> - The library name. The path to diversity file.

B<$file> - The path to diversity file.

I<Output Parameters>:

B<$output> - The DiversityOnlyIsotigs format output file.

=head2 B<generateIsotigList()> - Generate a list of isotigs that exist in each diversity file.

I<Input Parameters>:

B<$working_dir> - The working directory where the isotig_list file is located.

B<$file_list> - The isotig_list files to parse for isotigs.

I<Output Parameters>:

B<$isotig_file> - The isotigs list output file.


=head2 B<generateTechFile()> - Creates template technology file for OTU extraction used by the convertTech sub routine.

I<Input Parameters>:

B<$isotig_list_file> - The isotigs list format file.

B<$scaleValue> - The number of median reads.

I<Output Parameters>:

B<$tech_file> - The path to the technology file


=head2 B<convertTech()> - Converts the template files into technology (*.tsv) files usable by gene_spring.

I<Input Parameters>:

B<$working_dir> - The working directory where the output files are going.

B<$tech_file> - The path to the technology file

B<$files2Convert> - hashref which maps to the files to convert

I<Output Parameters>:

B<$status> - Status of the convertTech check

B<$genespring_dir> - The genespring directory

=head2 B<convertTechCheck()> - Checks the individual "diversityTech" output files for errors. If success returns successful message Otherwise kills program and prints the file/line number of the error.

I<Input Parameters>:

B<$genespring_dir> - The genespring directory

B<$tech_file> - The path to the technology file

I<Output Parameters>:

B<$status> - Status of the convertTech check


=head2 B<generateFastaMEGAN()> - Creates a MEGAN FASTA input file for each library

I<Input Parameters>:

B<$working_dir> - The working directory to construct a MOTHUR directory where the output files are going.

B<$divfile> - The path to diversity file.

B<$CDHitClustFnaFile> - The path to the clstr.fna file.

=head2 B<generateBLASTX()> -  Creates a BLASTx input file for the generateBlastxMEGAN sub-routine.   

I<Input Parameters>:

B<$assembly_Dir> - The assembly directory

B<$CDHitClustFnaFile> - The path to the clstr.fna file.

I<Output Parameters>:

B<$CDHitClustFnaBlastxFile> - The path to the clstr.fna.blastx file.


=head2 B<generateBlastxMEGAN()> - Creates a MEGAN BLASTX input file for each library

I<Input Parameters>:

B<$working_dir> - The working directory to construct a MOTHUR directory where the output files are going.

B<$divfile> - The path to diversity file.

B<$CDHitClustFnaBlastxFile> - The path to the clstr.fna.blastx file.



=head2 B<generateInputMOTHUR()> - Creates a MEGAN BLASTX input file for each library

I<Input Parameters>:

B<$working_dir> - The working directory to construct a MOTHUR directory where the output files are going.

B<$genespring_dir> - The genespring directory

=head2 B<processInputMOTHUR()> - Processes MOTHUR input files that consist of each library performs calc=nseqs-coverage-npshannon-simpson-sobs-chao and rarefaction data.

I<Input Parameters>:

B<$mothur_file> - The working directory to construct a MOTHUR directory where the output files is going.

B<$mothur_file_dir> - The genespring directory

=head2 B<generateInputUNIFRAC()> - Creates a UNIFRAC input file consisting of each library

I<Input Parameters>:

B<$working_dir> - The working directory to construct a UNIFRAC directory where the output file is going.

B<$genespring_dir> - The genespring directory

=cut


1;
