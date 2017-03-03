

Please note that this file is now deprecated. There has been a major refactoring of mPUMA to function through the use of Make for program execution. It is my intention to have corresponding documentation uploaded over the next couple of months. 

Matthew Links - 3March2017


README for mPUMA 

Dependencies that must be installed and configured prior to installing mPUMA. 

- PERL installation varies by version of Linux

- SFFTools and "GS De Novo Assembler" are distributed by Roche/454. To obtain a copy of them go to 454.com -> Products -> Analysis Software and request a free download of the software

- Trinity (optional) http://trinityrnaseq.sourceforge.net/
	- If you are looking to increase the kmer size from the default of 25 then you will need to edit Trinity.pl so set $IWORM_KMER_SIZE and $MIN_IWORM_LEN to be the new value for kmer size.
 
- Bowtie 2 http://bowtie-bio.sourceforge.net/bowtie2

- NCBI's standalone BLAST http://blast.ncbi.nlm.nih.gov

- SeqClean, cdbfasta and cdbyank http://compbio.dfci.harvard.edu/tgi/software/

- CDHit http://cd-hit.org

- Mothur http://www.mothur.org

- T-Coffee http://tcoffee.org

- FastTreeMP http://meta.microbesonline.org/fasttree/ install the multi-threaded 
version

- RDP classifier	http://sourceforge.net/projects/rdp-classifier/files/

The installation of these tools is beyond the scope of this document. 

A quick walkthrough of using mPUMA with example data
----------------------------------------------------

1. Create some directory to store projects you intend to access using SVN.

	mkdir ~/svnroot

2. Go into your svnroot

	cd ~/svnroot

3. Get a version of the mPUMA code.

You can checkout the current version of mPUMA from the souce repository.

	svn co https://svn.code.sf.net/p/mpuma/code mpuma

N.B. this will take a various amount of time depending on your Internet connection

4. Move into the mPUMA directory.

	cd mpuma

5. Figure out what your current directory is as you will need to set the mPUMA environment variable to this location.

	pwd

6. Set the mPUMA environment variable to the current directory. From now on $mPUMA refers to this location.

For BASH you would likely want to add the following to your ~/.bashrc

	export mPUMA=~/svnroot/mpuma

### Please note that from this point onwards the build process has changed to reflect the modified procedure within mPUMA that uses make to build all the dependancies for your analysis (Aug 2014)

7. Check default settings in the mPUMA Makefile. The configured paths, parameters, defaults are housed in $mPUMA/make/Makefile

	vi $mPUMA/make/Makefile

You should look at the variables that are being set with a default value by looking for ?= 

	Example: mPUMA ?= /usr/local/mpuma

This example shows that the mPUMA variable is being set to the path /usr/local/mpuma. This should be the location that you have pulled mPUMA down to via the svn co. 

You should pay special attention to the paths to each of the executables that mPUMA uses. If you aren't specifying a fully qualified path name then the assumption is that the executable is in your current path 

	# path to cdbfasta
	CDBFASTA = cdbfasta
	# path to cdbyank
	CDBYANK = cdbyank
	# path to CDHit-EST
	CDHITEST = /usr/local/cd-hit/cd-hit-est
	# path to CDHit for amino acid clustering
	CDHIT = /usr/local/cd-hit/cd-hit
	# path to blastall - needed for orientating OTUs
	BLASTALL = /usr/local/blast/bin/blastall
	# path to tCoffee
	TCOFFEE = t_coffee
	# path to FastTreeMP
	FASTTREE = FastTreeMP
	# path to bowtie2
	BOWTIE2 = /usr/local/bowtie2/bowtie2
	# path to bowtie2-build
	BOWTIE2BUILD = /usr/local/bowtie2/bowtie2-build

	# Location of the seqclean executable
	# N.B. SeqClean as a pre-compiled binary can usually only support up to 16 CPUs
	# the issue isn't acually in the PERL script SeqClean but rather in the psx command
	# it calls. It is recommended that for installations with > 16 CPUs you get the source
	# for psx and recompile it after you edit psx.h and change the line
	# define MAX_PROCS 16
	# to the desired number. Make sure that within the PERL SeqClean script it references
	# the correct psx executable. Once it does you should be able to pass through the CPU #
	SEQCLEAN = /usr/local/seqclean-19June2014/seqclean
	
As well as some presets for # of CPUs to use etc.

	NUM_CPU = 40
	# the non-redundant FASTA of reference sequences - used for simple mathing and C3
	CPNDB_NR = $(mPUMA)/reference_fasta/cpndb_nr_20140327
	CPNDB_NR_PEP = $(mPUMA)/reference_fasta/cpndb_nr_pep_20140327
	PRIMERS = $(mPUMA)/reference_fasta/primers-enumerated-3primeN.msf

8. Download the example data. N.B. this is ~ 200 MB in size so the time will be dependent on your Internet connection. 

one of the following should work

	curl -L  https://downloads.sourceforge.net/project/mpuma/Chaban_Links_Hill-Dog-Fecal-Microbiome.tar.gz > Chaban_Links_Hill-Dog-Fecal-Microbiome.tar.gz
	wget https://downloads.sourceforge.net/project/mpuma/Chaban_Links_Hill-Dog-Fecal-Microbiome.tar.gz

9. Uncompress and unpack the archive. 

	tar -zxvf Chaban_Links_Hill-Dog-Fecal-Microbiome.tar.gz

10. Create a library.spec file in the current working directory which specifies the location of the unpacked SFF files. 

	vi ~/mpuma/Canine-fecal-microbiome/Chaban_Links_Hill-Dog-Fecal-Microbiome/library.spec

The file should look something like:

DDS18   /path/GJBSLIE13-MID4.sff
DDS19   /path/GJBSLIE14-MID6.sff
DDS2    /path/GJBSLIE12-MID11.sff
DDS30   /path/GJBSLIE14-MID8.sff
DDS32   /path/GJBSLIE14-MID11.sff
DDS34   /path/GJBSLIE15-MID2.sff
DDS51   /path/GJBSLIE15-MID4.sff
DDS64   /path/GMH47UN08-MID6.sff
HDS18D  /path/GJBSLIE11-MID4.sff
HDS19   /path/GJBSLIE12-MID6.sff
HDS1fall        /path/GJBSLIE09-MID2.sff
HDS1spring      /path/GJBSLIE09-MID4.sff
HDS2fall        /path/GJBSLIE10-MID6.sff
HDS2spring      /path/GJBSLIE10-MID8.sff
HDS30   /path/GJBSLIE12-MID8.sff
HDS8    /path/GJBSLIE10-MID11.sff

The first column is the name you want the library to have. The second column is the path to the file. It is strongly suggested that you use full paths to these files. You can when working with other data have multiple SFF files listed per library.

==========HINT==========
Using a copy of the library.spec file, which was in the archived example data you might use vi to change the "path" in each row using the vi syntax similar to 

:0,$ s/\/path/\/home\/USERNAME\/mpuma\/Canine-fecal-microbiome\/Chaban_Links_Hill-Dog-Fecal-Microbiome/

Note - you CANNOT use a reference to your home directory which uses the tilde character "~". For some reason this is interpreted differently by open2 in Perl and leads to errors about now being able to find the specified file. 

========================

11. Using make to build a specified target. There are a couple of targets that you could choose to build: assemble, condense, abundance, biom, phylogeny. 

	assemble	This is the mPUMA de-novo assembly process. 
	condense	Depends on assemble and procedes with primer trimmining, clustering at 100% identity, and chimera checking (optional). 
	abundance	Depends on consense and will generate the table-of-contents file (total.toc) that allows the caluclation of OTU abundance
	biom		Depends on abundance and will produce a biom file (of type dense version 1)

Additionally there are a series of cleanup targets that can be used 

	clean_biom
	clean_abundance
	clean_condense
	clean_assembly 

Each of these aim to return to the state prior to the respective build target. 

Likely what you will want to do at this point is check to see if everything is working. So you set your sites on a biom file as the product

make -f $mPUMA/make/Makefile biom LIBRARY_SPEC=library.spec ASSEMBLY_METHOD=Trinity MAPPING_METHOD=bowtie2 ASSEMBLY_NAME=Trinity-bowtie2

C'est ca

Additional documentation - please be aware that the remaining comments are from the version of mPUMA that preceeded the make build process. These comments may or maynot apply to the current (Aug2014) version that is based on make and so please be conscious of that as you continue reading. Further updates to mPUMA are coming and will result in additional updates to this README too. 

Outputs
--------

There are 3 general partitions to the outputs of mPUMA 
-------------------------------------------------------------------

1. Assembly outputs contain the sequence outputs associated with deriving the OTUs (e.g. ~/mpuma/Canine-fecal-microbiome/31Oct2012)

By default mPUMA will create a datestamp in the format DDMMMYYYY to use for the assembly directory it creates (e.g. 31Oct2012). This directory is created as part of the calls to the 454/Roche assembly tools and you are likely interested in the files within the assembly subfolder (e.g. ~/mpuma/Canine-fecal-microbiome/31Oct2012/assembly). A description of all of the 454 outputs is beyond the scope of this document but you are encouraged to read the 454/Roche literature on gsAssembler (particularly the de novo assembler in cDNA mode) as it's described in the documentation you should have received from Roche/454 when you downloaded their tools. 

mPUMA does compile all of the input reads together (total.fna) and does call the classifier (total.fna.classified) on them was well as perform comparison to the reference database via wateredBLAST (total.fna.wateredBLAST.cpndb_nr_20110622). These analyses are not likely what most researchers will be interested in but they are used in assessment of assembly error (e.g. Specificity and Sensitivity). In most cases a researcher is interested in using the assembled OTUs and their abundance to evaluate hypotheses. 

Most of the files mPUMA creates are named in a highly linear fashion. The output of each mPUMA operation is usually responsible for the terminal suffix on the output filename. 

454Isotigs.fna - The initial set of OTUs which mPUMA creates are the direct output of the de novo assembly 

Seq-cleaning of the initial OTUs produces a report (454Isotigs.fna.cln), a log file (seqcl_454Isotigs.fna.log), and a set of OTUs which have had primers removed (454Isotigs.fna.seqclean). The seqcleaning step does use a fully enumerated sequence database of the primers (2.4M sequences for the degenerate universal cpn60 primers). In some cases there will be the odd primer sequence which SeqClean cannot detect and so there may be a token OTU sequence which still has the amplification primer present. 

Chimera checking via Chaban's Chimera Checker (C3). C3 is named for a chimera cleanup procedure proposed by Dr. Bonnie Chaban (Hill lab - UofS). First, the ends (150bps) from each of the seqclean'd OTUs are extracted (454Isotigs.fna.seqclean.c2e). The end sequences are then compared to a reference database (454Isotigs.fna.seqclean.c2e.wateredBLAST) and a report of the findings is produced (454Isotigs.fna.seqclean.c3report). If the user has specified that C3 chimera checking should be used (default behavior) then the end sequence matches will be examined to see if both ends match the same reference sequence and in the correct orientation. Any end-sequences that don't match correctly and in agreement with its matched end will cause an OTU to be flagged as possibly chimeric (454Isotigs.fna.seqclean.c3report.chimericIDs). Studies of ecological niches where you have good representation of sequences in the reference database should benefit greatly from this procedure. 

However, caution should be taken when the researcher knows that the reference database may be deficient in taxa of interest that will likely be detected in the experimental data. In such cases there will be an occurrence of false positive chimera calls where either end of an OTU may match different sequences in the reference database. Commonly each end will have a match to the same Genus but to different species with low percent identities. If you are working in a situation where you want to avoid using the C3 chimera calls this can be turned off by adding "-c 0" to the mPUMA pipeline and / or analysis scripts. If C3 chimera checking is off the analysis will still be performed but the chimera calls won't be used to exclude OTUs from your analysis. You would still have one level of chimera checking performed by gsAssembler when the OTUs were originally formed. 

454Isotigs.fna.seqclean.possibleChimerasRemoved is the result of the assembly + seqcleaning + chimera analysis (C3). Even if you have disabled C3 chimera calls this file will exist and is the next logical step in the mPUMA pipeline.

CD-Hit clustering is use twice in the mPUMA pipeline: once to collapse exactly identical sequences in the nucleotide space and once to collapse synonymous differences in the amino acid sequences of the OTUs. Both times CD-Hit is called by mPUMA, it performs clustering at 100% identity. 

In each case there is a clustering output which has a filename ending in .clstr (e.g. 454Isotigs.fna.seqclean.possibleChimerasRemoved.cd-hit.clstr or 454Isotigs.fna.seqclean.possibleChimerasRemoved.oriented.aa.cd-hit.clstr). That clustering output is parsed to create a clustering report which has a name ending in .report (e.g. 454Isotigs.fna.seqclean.possibleChimerasRemoved.cd-hit.clstr.report or 454Isotigs.fna.seqclean.possibleChimerasRemoved.oriented.aa.cd-hit.clstr.report). Each clustering report is a tab delimited text file which contains the columns OTU name of the longest member of the cluster, Length of the longest member of the cluster, # of member OTUs which can be safely collapsed into this larger sequence, a comma separated list of the shorter OTUs which have been collapsed. The clustering report is crucial for tracking each experimental read into its assembly OTU and then into the representative OTU for that cluster. 

A final FASTA is produced for each round of clustering which has a filename ending .cd-hit.clstr.fna (e.g. 454Isotigs.fna.seqclean.possibleChimerasRemoved.cd-hit.clstr.fna  or 454Isotigs.fna.seqclean.possibleChimerasRemoved.oriented.aa.cd-hit.clstr.fna). 

When mPUMA completes the clustering of the OTUs in the nucleotide sequence space it carries out a number of analyses on the non-redundant set of OTUs. Classification using a Bayseian Classifier (from RDP codebase but trained on cpn60 sequences) produces output (454Isotigs.fna.seqclean.possibleChimerasRemoved.cd-hit.clstr.fna.classified) which can be used to create percentage abundance composition profiles (CLASSIFIER_PROFILES) or for input to MEGAN (FILES_FOR_MEGAN). 

Searches against a reference database are done via wateredBLAST (454Isotigs.fna.seqclean.possibleChimerasRemoved.cd-hit.clstr.fna.wateredBLAST.cpndb_nr_20110622) to provide a simple nearest neighbour match, which could be used in downstream analyses as a label. 

Comparison to a reference database via BLASTX produces output (454Isotigs.fna.seqclean.possibleChimerasRemoved.blastx) that can be used to identify a likely frame for protein translation. In the case of protein coding targets (protein coding steps can be turned off with "-protein 0" in calls to mPUMA-pipeline-analysis.pl) this likely frame is used to orient (454Isotigs.fna.seqclean.possibleChimerasRemoved.oriented) the nucleotide sequences and then translate them (454Isotigs.fna.seqclean.possibleChimerasRemoved.oriented.aa). 

N.B. the phylogenetic analysis of the nucleotide sequences requires that the OTU sequences are oriented so even if you have turned off the protein coding analyses you will still need to have a homology search performed to determine the sequence orientation prior to aligning the non-redundant nucleotide sequence OTUs. 

mPUMA accomplishes an automated phylogenetic analysis for both the nucleotide and amino acid OTUs. In each case the multiple sequence alignment is performed using tCoffee and the phylogenetic analysis is performed via FastTree. The outputs from tCoffee will include the file names ending in .dnd and .aln. Where as the output from FastTree (the actual phylogenetic tree) will end in .FastTree. 

When mPUMA uses gsAssembler for read to OTU tracking (N.B. this is NOT the default anymore).
The way mPUMA calls gsAssembler ensures that there is an "ace" directory created which contains 1 ace file per assembled OTU. This is the largest output commonly produced in running mPUMA. Once mPUMA has processed each of the ace files created in the de novo assembly step they will be compressed via gzip to conserve space (e.g. isotig00001.ace.gz). In processing the OTU files mPUMA figures out which reads are in each ace file (isotig00001.ace.reads) and then traces all reads together into a table of contents file for the assembled OTUs (ace.toc). The table of contents file for the assembly is a tab delimited text file which contains 2 columns, OTU and read.

N.B. Since gsAssembler is run in cDNA mode by mPUMA there is technically the chance for a read to be used in the assembly of 2 OTUs. To gsAssembler this is a case where 1 read fits into isoforms for a gene. The conceptual parallel to OTU assembly is a case where the assembler is building a partial (e.g. 5'-end) OTU sequence and a full length one and has a read that can be used (correctly) to build both. The way that this can cause confusion for the researcher (and it seems to happen to all of us the first time do this) is when you try to tally up exactly the number of reads produced and naively assume that each read can go to only 1 place in the assembly (which is the false assumption we all make). To track exactly all read fates is complicated even in this type of assembly (overlap-layout-consense) as many algorithms (gsAssembler included) trash some sequences for a myriad of QC reasons even before they start aligning them. 

When mPUMA uses Bowtie 2 for read to OTU tracking (this is the default behavior as of March 2013), mPUMA employs Bowtie 2 to map all of the Reads onto the assembled OTU sequences. Prior to actually doing the mapping there will be a conversion of any of the SFFs as necessary to Fastq. The inputs for the Bowtie 2 exercise will reside within a folder called fastq. The outputs from the Bowtie 2 mapping for each of the input files will be placed into an output folder within the assembly directory called sam. The sam files are parsed to create toc file sam.toc which traces each read to the assembled OTUs and then these are compiled with information on the cleanup steps from OTU formation (seqcleaning, and cd-hit) to record the final resting place for each read in the total.toc.

Option - Running Trinity as the assembler

If you chose to use trinity as the assembler then the output from the assembly will be called Trinity.fasta (instead of 454Isotigs.fna from gsAssembler). Then all of the seq-cleaning and cd-hit steps will add suffixes to the file name analogous to the ones described above for when gsAssembler is used to assemble the OTU. 


In order to facilitate read tracking through the entire assembly process, mPUMA analyzes the results of de novo assembly, seqcleaning, chimera checking, and clustering. A final map, which indicates where each read ends up in the mPUMA output, is produced as a master table of contents (total.toc). The output is a comma separated text file that has the columns readID, final resting place of the read, rationale. This makes it straightforward to do some interactive QC'ing of the process. For example suppose you wanted to assess whether the chimeric calls by C3 were affecting a large proportion of your reads you could to something like 

cat ~/mpuma/Canine-fecal-microbiome/31Oct2012/assembly/total.toc | cut -f 3 -d ',' | grep -v Rationale | sort | uniq -c | sort -nrk 1

which will break down for you (by # of reads) what the fate of the reads were in the mPUMA analysis. Ideally you would want to see all of your data in "CDhit-OTU-was-redundant" or "CDhit-OTU-was-non-redundant". In practice however you should expect that some proportion of your data lost by QC assessment in the assembly step (gsAssembler-left-out) or in the SeqCleaning of the sequences (seqClean-trashed). If you are seeing a high (sorry, but this is a subjective measure at the moment) proportion of your data being flagged as being in "Chimeric-OTU" then you should consider not using C3 to exclude those data, and make sure you convince yourself that any interesting OTUs which you want to carry on investigating are not chimeric by additional methods (e.g. manual curation, bellerophon, etc.). 

==========HINT==========
If you have run the entire mPUMA pipeline once and then realized you want to exclude the C3 chimera calls, or make other changes to the way the post-assembly scripts were run, you can do the following. Note that the assembly will NOT be re-computed, only the post-assembly analysis.

a) Decide whether you want to keep a backup copy of the outputs from the original analysis. If you want to compare results of different post-assembly analysis procedures then you should probably copy the nucleotide and protein output folders to another location on disk. 

b) Delete all the mPUMA files from C3 onward. WARNING - make sure you have the period in the command ".". If you left it out you would loose your SeqCleaning results too and that can be painful to have to recalculate. 

	rm ~/mpuma/Canine-fecal-microbiome/31Oct2012/assembly/454Isotigs.fna.seqclean.*

or if Trinity was the assembly approach

	rm ~/mpuma/Canine-fecal-microbiome/31Oct2012/assembly/Trinity.fasta.seqclean.*

c) Delete the total.toc table of contents. 

	rm ~/mpuma/Canine-fecal-microbiome/31Oct2012/assembly/total.toc

d) Remove the nucleotide and protein output directories. LAST WARNING - if you wanted to compare between analyses then make sure you've archived the existing versions.

	rm -rf ~/mpuma/Canine-fecal-microbiome/nucleotide
	rm -rf ~/mpuma/Canine-fecal-microbiome/protein


e) Re-run the post-assembly analyses while avoiding the use of the C3 chimera calls (the report is still generated, just not used) 

	$mPUMA/bin/mPUMA-pipeline-analysis.pl -c 0 -w ~/mpuma/Canine-fecal-microbiome -l ~/mpuma/Canine-fecal-microbiome/library.spec -n 31Oct2012 > ~/mpuma/Canine-fecal-microbiome/log.noC3 2>&1

========================

Outputs for nucleotide and protein sequence OTUs
------------------------------------------------

2. Outputs for use with other tools based on the nucleotide sequence of the OTUs (e.g. ~/mpuma/Canine-fecal-microbiome/nucleotide)
&
3. Outputs for use with other tools based on the protein sequence of the OTUs (e.g. ~/mpuma/Canine-fecal-microbiome/protein)

For the outputs associated with either the nucleotide or protein sequence of the OTUs there are files broken down as follows:
	
libraries
		
A directory containing a sub-directory for each library named in library.spec

For each library there will be:
	reads.ids (A text file containing the list of every experimental sequence associated with that library.)
	subsample.ids (A text file containing a list of sequences from the library chosen at random. The total number will be either user defined or the minimum library size.)
	diversity (A text file containing abundance information on the OTUs detected in this library. Abundance of each OTU is determined from the sequence IDs in subsample.ids tracked through the assembly procedure. The file is tab delimited and contains columns OTU, abundance (percentage), scaled abundance (scaled to 10000 reads), absolute abundance, Best match by wateredBLAST. The 5th column which contains a string for the best match has a comma separated value of the ID + description of the match, % identity of the alignment, length of the alignment, orientation.

CLASSIFIER_PROFILES
A directory containing a text file for each library named LibraryName.taxonomicLevel. Each of the text files describes the library in terms of its taxonomic composition (100%) to the particular taxonomic level. The default taxonomic level is phylum and is controlled by the -t parameter to $mPUMA/bin/mPUMA-pipeline-analysis.pl. The profiles are based on the Bayesian classifier results for the OTUs and weighted according to the abundance coming from the reads in the libraries subsample.ids file. If you are interested in the original classifier outputs you should look within the assembly directory.

FILES_FOR_MOTHUR
A directory containing the input for mothur (freq_mothur.txt), the rabund outputs, and the outputs in a subfolder called MOTHUR_OUTPUT_FILES
		- MOTHUR_OUTPUT_FILES contains a rarefaction output for each library and the overall ecological statistics in freq_mothur.groups.summary (e.g. ~/mpuma/Canine-fecal-microbiome/nucleotide/FILES_FOR_MOTHUR/MOTHUR_OUTPUT_FILES/freq_mothur.groups.summary)

FILES_FOR_MEGAN
A directory containing a file for each library which can be imported into MEGAN as standalone RDP formatted data. Each of these files is a weighted interpretation of the OTU classifier data where each libraries subsample.ids file is used to infer the abundance of each OTU. For each OTU its classifier result is repeated in the input for MEGAN n times where n was the abundance in that library. The first column of the RDP format contains the ID of the sequence being classified. To track abundance the id is given a suffix of "_n" where n is the nth time that OTU was observed in that library. There are a number of reasons the input for MEGAN is formatted in this way: 1) the use of the RDP format is much smaller and more manageable then providing FASTA, or BLASTX results to MEGAN and 2) when we have experimented with providing input to MEGAN based on processed BLASTX results (i.e. 3 column format) there has been no way to trace the assignments made internally by MEGAN.

FILES_FOR_UNIFRAC
A directory containing 3 files:
	Categories.txt (This is a rough sketch of what your metadata for unifrac will need to look like. The library names are pulled from library.spec. http://bmf2.colorado.edu/fastunifrac/tutorial.psp)
	FastTree (This is the phylogenetic tree which was generated from tCoffee and FastTree for these data. Within this directory this is a symbolic link to the actually file which resides in the assembly directory.
	sampleIDmappingFile.txt (This is the sample abundance data which unifrac will use. http://bmf2.colorado.edu/fastunifrac/tutorial.psp)

To use these results with Unifrac the user MUST make some modifications. The Categories.txt file needs metadata entered into it for unifrac to be able to analyze your experiment. Uploading of these files to the Unifrac web resources typically requires that the user compress them with gzip (e.g. gzip Categories.txt). 
		
FILES_FOR_GENESPRING
A directory containing 1 text file for each library. Each file is a tab delimited text file with the columns OTU, Percent abundance, Abundance scaled to 10000 reads / library, Actual count, Label. These files can be used with analyses which compare abundance across a series of samples / libraries. This can be done using a tool like GeneSpring. When using GeneSpring you need to create a custom 1-colour technology using any of the files from this directory. By creating a custom technology you can get GeneSpring to interpret the abundance data from whichever column you feel is appropriate as if it were a fluorescence reading from a microarray study. These files can also be used with R to perform similar analyses. Look for separate mPUMA documentation on these types of analyses in the near future. 

OTU-list
This is simply a text file which lists all of the OTUs which were detected in the assembly. Technically this is taking into account the down-sampling exercise (i.e. it is based on the subsample.ids files for all libraries), so in practice this may actually exclude some extremely rare OTUs which were not sampled when downsampling your libraries at random. 

OTU-list.tech.tsv
This is a skeleton of what the files within the FILES_FOR_GENESPRING folder will look like. This is likely of no use unless you are a developer and this file may disappear from future releases. 

Advanced approaches & pointers
------------------------------------
Host genome cleanup

When you are studying a microbiome as it is linked to a specific host you may want / need to cleanup the data for host sequences. In most cases the host is likely a eukaryote and if there is a genome sequence (of any state) available you can use it to clean up your data prior to using mPUMA. The chief thing one would aim to clean up for would be off-target amplifications. When degenerate PCR primers are used it is commonplace to get some unintended amplifications from templates which are complementary to the degenerate primers but are not the intended barcode / gene. The suggested procedure for mPUMA involves masking all cpn60 based locations in the host's genome and then using the masked non-cpn60 host genome to screen data prior to using mPUMA. This does mean that we are suggesting that you keep the cpn60 regions of the host genome in the analysis. 

1. Find the locations of the cpn60 or cpn60 pseudogenes in the host genome. Performing a tblastn search with a cpn60 UT amino acid sequence from cpnDB easily accomplishes this. It's suggested that you find the most closely related organism in cpnDB and use tblastn to fish into the host genome and pull out the best match. This should easily identify for you an HSP, which corresponds to a cpn60 from the host genome. 

2. Extract the sequence of the host's cpn60 sequence. You can do this easily with extractseq (EMBOSS) to create a fasta file containing the sequence of a host cpn60 gene. 

3. Given that it is commonplace for higher eukaryotes (plants and mammals) to have multiple copies of cpn60 and even some pseudogenes we recommend that you identify these by performing a second BLAST search into the host genome but this time using the host sequence from the previous step as bait.

4. Create a text file containing the genomic locations (just use the 5' and 3' extremes of the HSPs from BLAST, you don't need to worry about splicing of introns). The file should be formatted with 3 columns, whitespace delimited where the 1st column is the sequence ID of the sequence containing a cpn60 gene or pseudogene. The second column contains the start of the cpn60 and the third column contains the end of the match. 

5. Using $mPUMA/bin/multi_masker.pl you should now be able to mask your host genome with a command similar to the one below. The masking_file is the whitespace delimited file from the previous step. 

$mPUMA/bin/multi_masker.pl -i host_genome.fasta -o host_genome_masked_cpn60.fasta -m masking_file

6. The actual screening vs. the host genome is performed on Fastq formatted input (HINT: if you need to convert SFF to Fastq check out $mPUMA/bin/sffToFastq.pl) using $mPUMA/bin/screen_host_genome_fastq.pl. Internally bowtie2 is called so you will need to build a bowtie2 database of the masked genme (HINT: see bowtie2-build command). The screening can be done on Fastq files which are compressed via gzip or the uncompressed versions and would use a command similar to 

$mPUMA/bin/screen_host_genome_fastq.pl -h host_genome_bowtie2_base -c numCPUsAllowed -o /output/Directory/HOSTSCREENED/ -z 0 file0 file1,? fileN

this will create files in /output/Directory/ which have the EXACT same basename as the files you specified for screening.

N.B. IT IS CRUCIAL THAT YOU NOTE THAT THE BASENAMES OF THE FILES BEFORE AND AFTER SCREENING WILL BE IDENTICAL. IT IS RECOMMENDED THAT THE OUTPUT DIRECTORY CONTAIN A REFERENCE IN ITS NAME TO THE FACT THAT THIS ARE MODIFIED FILES!!!

There are a couple of host genomes which have already been masked for the cpn60 locations available on the sourceforge website. At the moment these included human, tomato and pig. We would suggest you place these datasets in a place like $mPUMA/reference_fasta/. If you have additional host genomes which you have masked for cpn60 please let us know if you are interested in depositing those for others to use and we can try and get them into the sourceforge repository. 

PLEASE BE AWARE THAT THE USE OF A HOST GENOME FOR SCREENING OF OFF-TARGET AMPLIFICATIONS LIKELY REQUIRES AN ACKNOWLEDGEMENT TO THE GROUP WHO RELEASED THE SEQUENCE.

Normalizing proportional abundance to a specific level

Suppose that you have carried out an estimation of the # of organisms in each of your samples independently of your sequencing experiment and you wish to scale the proportional abundances to this level. E.g. you calculated 16S rRNA copies / g and are using this as an estimate of microbial load.

N.B. THIS IS FUNDAMENTALLY DIFFERENT FROM RANDOM SAMPLING TO THE SMALLEST LIBRARY SIZE. SUCH RANDOM SAMPLING IS SUGGESTED FOR THE CALCULATION OF ECOLOGICAL PARAMETERS SUCH AS ALPHA AND BETA DIVERSITY ESTIMATES AS WELL AS EVENNESS. THE PROCEDURE DISCRIBED HERE IS INTENDED FOR CASES OF PATTERN DISCOVERY (E.G. CLUSTERING).

1. Create a new directory to store these data. 
2. For each library you wish to scale you should be able to use a command like the following to scale the reads from LIBRARY-A to 10M.

$mPUMA/bin/generate_diversity_from_toc.pl -t /path/to/assembly/total.toc -i /path/to/nucleotide/libraries/LIBRARY-A/reads.ids -o /path/to/new/directory/LIBRARY-A-scaled-to-10M.diversity -c /path/to/assembly/454Isotigs.fna.seqclean.cd-hit.clstr.report -scale 10000000


Creating a list of OTUs based on a prevalence threshold.

You can use the following to print a list of OTUs which have a prevalence of at least INTEGER (# of libraries which it has to occur in, at least this many). 

$mPUMA/bin/find-prevalent-OTUs.pl -o 0 -h 0 -p 7 File0 file1 file2 ? fileN > min-prevalence-7-of-11-libraries.ids 

The options for -o specifices a 0-index column within the files (whitespace delimited) which contains the OTU label, -h is a boolean for whether to discard the first line as a header and -p is the Interger value of the minimum prevalence. 


Working with diversity data from mPUMA in R

If you want to use R for various downstream analyses you probably will want to use $mPUMA/bin/generate_analysis_files_R.pl to generate you a matrix formatted file which you can use within R. So assuming you have created a directory which contains your diversity files (filenames ending in .total-diversity, perhaps created from a scaling approach as above) you could use a command like the following to parse all of the diversity files in /path/to/new/directory/ and create a suitable input file for R.

$mPUMA/bin/generate_analysis_files_R.pl -i min-prevalence-7-of-11-libraries.ids -d /path/to/new/directory/ -s .total-diversity -c 2 -o R/input-for-R-minPrevalence7-scaled-10M-.txt


Known gotchas
-------------

PROBLEM - gsAssembler errors. We have noticed that gsAssembler can generate errors which seem to arise from some type of race condition. These are most commonly encountered as problems where mPUMA fails to run the post assembly analyses. Investigating the last couple of lines of an mPUMA run log shows the following error message

Error:  Initialization failure in ThreadMutex constructor.        
        Please report this error to your customer support representative.

Error:  An internal error (assertion failure) has occurred in the computation.
        Please report this error to your customer support representative.

solution - rerun the analysis. Delete the assembly directory and rerun the pipeline. As this seems to be a random or race condition it does not tend to reoccur. I know that this is not reassuring when using a computer but it is a feature / nuance of gsAssembler and may be version specific.


PROBLEM - some error occurs and you want to re-run the pipeline but not re-assemble everything. If you don't specify the assembly name (-n) when running the $mPUMA/bin/mPUMA-pipeline.pl command it will by default use a datestamp in the format DDMMMYYYY (e.g. 31Oct2012). 

Solution - If you had some kind of error occur post assembly, and at least 1 day has passed since the assembly ran then you likely need to add a -n parameter (e.g. "-n 31Oct2012") to the $mPUMA/bin/mPUMA-pipeline.pl command line. Alternatively if you know that the assembly is ok then you probably should just call the post assembly analysis portion of the pipeline directly $mPUMA/bin/mPUMA-pipeline-analysis.pl.


PROBLEM - you either know or suspect that your ecological niche is not well represented in cpnDB. If this is the case and you use the putative chimera calls from C3 you could easily end up excluding huge proportions of your data as false chimeras. 

Solution - You can turn off the use of the C3 chimera calls on the command line via -c 0. If you do this the first time you call the pipeline everything will be fine. C3 will be called and you will get the reports generated about chimeras they just won't be used to explicitly exclude OTUs from subsequent steps. If you have already run the assembly once and need to
re-run the post assembly steps then there are a couple of cleanup steps you need to take in your assembly directory before you re-run the analysis. Failure to do this may leave the C3 chimera calls in effect
              
	- delete assembly/total.toc
	- delete all of the files post seqcleaning assembly/454Isotigs.fna.seqclean.*
	- delete all of the outputs. From within both the nucleotide and
	  protein output folders delete
		CLASSIFIER_PROFILES
		FILES_FOR_GENESPRING
		FILES_FOR_MEGAN
		FILES_FOR_MOTHUR
		FILES_FOR_UNIFRAC
		OTU-list
		OTU-list.tech.tsv


PROBLEM - you want to change the 1st pass OTU file from 454Isotigs.fna to 454AllContigs.fna

Solution - You can specify the new file using -i and re-running the post-analysis steps. You MUST remove some of the files prior to re-running the post-analysis 
	- delete assembly/total.toc
	- delete assembly/ace.toc
	- delete all of the outputs. From within both the nucleotide and protein output folders delete
		CLASSIFIER_PROFILES
		FILES_FOR_GENESPRING
		FILES_FOR_MEGAN
		FILES_FOR_MOTHUR
		FILES_FOR_UNIFRAC
		OTU-list
		OTU-list.tech.tsv


PROBLEM - you actually want to use ALL reads and so you want to NOT down-sample the abundance. 

Solution - Set the number specified with "-downSample #" to be larger than
your largest library. That will cause the subsampling routine to select every read for inclusion in the calcuations. If you have already run the post assembly analysis then you will need to remove all of the outputs. So remove both the nucleotide and protein output directories entirely.


PROBLEM - you are getting permission errors. 

Solution - check that you can execute each of the 3rd party tools.
Configuration of those tools is out of scope for the mPUMA documentation. If you are having problems with something to do with reading / writing files check permissions on the directories and files. Also check and ensure that the working directory parameter (-w) is a fully qualified path. Relative paths may cause unexpected behaviour.
