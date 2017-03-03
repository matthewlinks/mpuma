all: abundance

# .SECONDARY side-effect: if target T depends on D and D is deleted, D won't be remade by `make T`.
# see https://www.gnu.org/software/make/manual/html_node/Special-Targets.html and 
# https://www.gnu.org/software/make/manual/html_node/Chained-Rules.html#Chained-Rules
.SECONDARY:

# assumption is the environment variable mPUMA is defined to point to the mPUMA installation 
mPUMA ?= /usr/local/mpuma

# some defualts for building .biom files
GENENAME ?= cpn60
ISPROTEIN ?= 1
PRETTYPRINT ?= 0

# Where is mPUMA going to produce outputs
CWD = $(shell pwd)
#These variables should be over-rideable via command line
# Where is most of the work going to happen
WORKING_DIR ?= $(CWD)

# Where is the configuration file?
LIBRARY_SPEC ?= $(CWD)/library.spec.sff

# Which assembly method? gsAssembler | Trinity
# Recomendation is for Trinity. gsAssembler is only appropriate if you have all SFF input.
# gsAssembler does seem to perform slightly better than Trinity for a synthetic community
# but it does not scale as well (Overlap-Layout-Consense vs. DeBruijn-Graph)
ASSEMBLY_METHOD ?= Trinity
# for gsAssembler
MIN_LEN ?= 137
MIN_ID ?= 90
# for Trinity the k-mer is determined in Trinity
# Recomended edits for Trinity are to set the Inchworm kmer to 31 
# my $IWORM_KMER_SIZE = 31;
# my $MIN_IWORM_LEN = $IWORM_KMER_SIZE;
# This is currently the maximum kmer size supported. There have been rewrites to some 
# of Trinity's dependancies to support larger kmers but these don't include all the 
# packages and so we can still only go up to a 32-bit word for the kmer

# Where are the enumerated primers used in the PCR (for seqcleaning)
# if you have the IUPAC coded version of your primers you can use 
# $mPUMA/bin/enumerate.pl to expand these out fully
PRIMERS	?= $(mPUMA)/reference_fasta/primers-enumerated-3primeN.msf

# Location of the seqclean executable
# N.B. SeqClean as a pre-compiled binary can usually only support up to 16 CPUs
# the issue isn't acually in the PERL script SeqClean but rather in the psx command 
# it calls. It is recommended that for installations with > 16 CPUs you get the source
# for psx and recompile it after you edit psx.h and change the line 
# define MAX_PROCS 16
# to the desired number. Make sure that within the PERL SeqClean script it references 
# the correct psx executable. Once it does you should be able to pass through the CPU #
SEQCLEAN = seqclean

# Do you want to use Chaban's Chimera Checker? 
# Recomendation: use when you are confident in having a representative FASTA for this ecological niche
# improper use of C3 could result in you flagging real OTUs as chimeras
USEC3 ?= 0

# What method would you like to use for read tracking? gsAssembler | bowtie2
# Recomendation is for bowtie2
MAPPING_METHOD ?= bowtie2

# a directory that is going to hold the assembled results
# this will contain a folder with symlinks to the data and an actual assembly folder with the results of building things
ASSEMBLY_NAME ?= $(ASSEMBLY_METHOD)-$(USEC3)
# The thinking here is that the assembly method defines everything downstream. 
# It is possible to have both forms of read tracking when we consider gsAssembler
# and so its not worth trying to capture the mapping method in the directory name
# USEC3 is embedded in the directory name because the C3 steps are run regardless 
# the only thing that changes is whether the putative chimeras are removed prior
# to the rest of the analysis continuing

# additional defaults that can be overridden
# make noise?
VERBOSE = 1
# how many cpu cores can be used at once? Warning a vanilla SeqClean install is limited to 16 
# if you are doing > 16 see the SeqClean comments above
NUM_CPU = 40
# the non-redundant FASTA of reference sequences - used for simple mathing and C3
CPNDB_NR = $(mPUMA)/reference_fasta/cpndb_nr_20140327
CPNDB_NR_PEP = $(mPUMA)/reference_fasta/cpndb_nr_pep_20140327
# path to cdbfasta
CDBFASTA = cdbfasta
# path to cdbyank
CDBYANK = cdbyank
# path to CDHit-EST
CDHITEST = cd-hit-est
# path to CDHit for amino acid clustering
CDHIT = cd-hit
# path to blastall - needed for orientating OTUs
BLASTALL = blastall
# path to tCoffee
TCOFFEE = t_coffee
# path to FastTreeMP
FASTTREE = FastTreeMP
# path to bowtie2
BOWTIE2 = bowtie2
# path to bowtie2-build
BOWTIE2BUILD = bowtie2-build

# Winnowing is a process where OTU that are below a WINNOW_CUTTOFF of anything in CPNDB_NR
# will be exlcluded from the final set of OTU. This 55% identity cuttoff has come from 
# years of looking at alignments of cpn60UTs in both the DNA and amino acid sequence spaces 
WINNOW ?= 1
WINNOW_CUTOFF ?= 55

# the path to where most stuff lives below
# this unfortunately was termed the ASSEMBLY_DIR
# when it actually contains within it a subdir called assembly
# you can thank Roche/454 for the idea of an 
#    Assembly_project/assembly 
# structure. Oi vei!
ASSEMBLY_DIR = $(WORKING_DIR)/$(ASSEMBLY_NAME)
SFF_DIR = $(ASSEMBLY_DIR)/sff
FASTQ_DIR = $(ASSEMBLY_DIR)/fastq
SAM_DIR = $(ASSEMBLY_DIR)/assembly/sam

# Multi-sequence FASTA file that is the product of the assembly effort
# Trinity.fasta or 454Isotigs.fna 
# One alternative idea would be to use 454AllContigs.fna for gsAssembler...
ASSEMBLED_OTU =
ifeq ($(ASSEMBLY_METHOD),gsAssembler)
	ASSEMBLED_OTU = $(ASSEMBLY_DIR)/assembly/454Isotigs.fna
else ifeq ($(ASSEMBLY_METHOD),Trinity)
	ASSEMBLED_OTU = $(ASSEMBLY_DIR)/assembly/Trinity.fasta
else
$(error Acceptable values for ASSEMBLY_METHOD are gsAssembler|Trinity not [$(ASSEMBLY_METHOD)])
endif

# Logic for the C3 processing
# Need to conditionally use the chimera check based on whether USEC3 is defined
OTU_PRE_CDHIT =
ifeq ($(USEC3),1)
	# the user does want to use the chimera calls and so the input to CD-Hit-EST should be  
	OTU_PRE_CDHIT = $(ASSEMBLED_OTU).seqclean.possibleChimerasRemoved
else ifeq ($(USEC3),0)
	# the user does NOT want to use the chimera calls and so the input to CD-Hit-EST should be  
	# $(ASSEMBLED_OTU).seqclean
	OTU_PRE_CDHIT = $(ASSEMBLED_OTU).seqclean
else
$(error Acceptable values for USEC3 are 0|1 not [$(USEC3)])
endif

# Need to conditionally use winnowing based on value of WINNOW
FINAL_OTU_SET = 
ifeq ($(WINNOW),1)
	# the user does want to winnow the data and exclude things that are nothing like cpn60UTs
	FINAL_OTU_SET = $(OTU_PRE_CDHIT).cd-hit.winnowed
else ifeq ($(WINNOW),0)
	# the user does NOT want to winnow the results
	FINAL_OTU_SET = $(OTU_PRE_CDHIT).cd-hit
else
$(error Acceptable values for WINNOW are 0|1 not [$(WINNOW)])
endif

# LIBRARY_SPEC contains the information on the inputs
# the file has a simple 2 column format where the 1st column is a library name (watch compatibility with Mothur)
# the 2nd column is a path to the inputFile 
# inputFile can end in sff|fastq|fq|fastq.gz|fq.gz
# All inputs
INPUT_FILES = $(shell cat $(LIBRARY_SPEC) | cut -f2)
# only the .sff files
SOURCE_SFF_FILES = $(shell cat $(LIBRARY_SPEC) | cut -f2 | grep -e '.sff$$')
# the fastq data regardless of compression status
SOURCE_FQ_FILES = $(shell cat $(LIBRARY_SPEC) | cut -f2 | grep -E '\.(fastq|fq)(\.gz)?$$')

# This is a template that allows the SFF files to be compressed to Fastq
# the conversion is needed when bowtie2 is being used for read tracking 
# and by necessity also when Trinity is the assembly method
FASTQ_FILES = 
# convert SFF to fq
define make-fq-sff
$(FASTQ_DIR)/$(notdir $(1)).fq: $(1)
	mkdir -p $(FASTQ_DIR)
	$(mPUMA)/bin/sffToFastq.pl -i $$< -o $$@ -z 0
FASTQ_FILES += $(FASTQ_DIR)/$(notdir $(1)).fq
endef
$(foreach x,$(SOURCE_SFF_FILES),$(eval $(call make-fq-sff,$x)))

# In order to perform read tracking via bowtie2 we create a fastq directory that 
# just contains a series of symbolic links to the Fastq data (compressed or not)
# that was specified in LIBRARY_SPEC
# N.B. this is the same directory with the converted SFFs will be put 
# so the directory actually can contain symbolic links to Fastqs (that can be gziped)
# as well as actual fastq data that is the result of conversion from SFF
define make-fq-fq
$(FASTQ_DIR)/$(notdir $(1)): $(1)
	mkdir -p $(FASTQ_DIR)
	ln -s $$< $$@
FASTQ_FILES += $(FASTQ_DIR)/$(notdir $(1))
endef
$(foreach x,$(SOURCE_FQ_FILES),$(eval $(call make-fq-fq,$x)))

# When doing the read tracking via bowtie2 we need a list of all of the possible SAM files
ALL_SAM_FILES =
ifneq ($(SOURCE_SFF_FILES),)
	ALL_SAM_FILES += $(foreach x,$(notdir $(SOURCE_SFF_FILES)),$(SAM_DIR)/$x.fq.sam)
endif
ifneq ($(SOURCE_FQ_FILES),)
	ALL_SAM_FILES += $(foreach x,$(notdir $(SOURCE_FQ_FILES)),$(SAM_DIR)/$x.sam)
endif

# The selection of MAPPING_METHOD determines which file is used to track 
# Reads <=> OTU
READ_TOC = 
ifeq ($(MAPPING_METHOD),bowtie2)
	READ_TOC = $(ASSEMBLY_DIR)/assembly/sam.toc
else ifeq ($(MAPPING_METHOD),gsAssembler)
	READ_TOC = $(ASSEMBLY_DIR)/assembly/ace.toc
else
$(error Asseptable values for MAPPING_METHOD are gsAssembler|bowtie2 not [$(MAPPING_METHOD)])
endif

# build sam.toc from all of the sam files
$(ASSEMBLY_DIR)/assembly/sam.toc: $(ALL_SAM_FILES)
	$(mPUMA)/bin/create_ace_toc_SAM.pl -s 1 -o $@ $(SAM_DIR)/*.sam

# build ace.toc from all of the ace files that are produced by gsAssembler
# this depends on the assembly having finished
# 2 possible cases where this could break
# a) somehow runProject is called without acedir
# b) somehow ASSEMBLY_METHOD is Trinity an somewhere that ace.toc gets asked for
# both of those would be non-sensical 
$(ASSEMBLY_DIR)/assembly/ace.toc: $(ASSEMBLED_OTU)
	$(mPUMA)/bin/make_toc.pl -d $(ASSEMBLY_DIR)/assembly/ace -o $@ -n $(NUM_CPU) 

# Dump all of the input read data to a file
$(ASSEMBLY_DIR)/assembly/single.fa: $(INPUT_FILES)
	mkdir -p $(ASSEMBLY_DIR)/assembly
	$(mPUMA)/bin/dump_fasta_from_sff_or_fastq.pl $(INPUT_FILES) -o $@

# There is a vestage between having total.fna (mPUMA thing) and single.fa (Trinity thing)
# so just use a symbolic link to support both
$(ASSEMBLY_DIR)/assembly/total.fna: $(ASSEMBLY_DIR)/assembly/single.fa
	ln -s $(ASSEMBLY_DIR)/assembly/single.fa $(ASSEMBLY_DIR)/assembly/total.fna

# need to concatenate all of the fq and fq.gzs into a file
# just use cat and zcat
$(ASSEMBLY_DIR)/assembly/total.fq: $(FASTQ_FILES)
	mkdir -p $(ASSEMBLY_DIR)/assembly
	$(mPUMA)/bin/dump_concatenated_fq_or_fqgz.pl $(FASTQ_FILES) -o $@

# Run Trinity on the basis of the total input being in place
# See note above about kmer size in Trinity.pl 
$(ASSEMBLY_DIR)/assembly/Trinity.fasta: $(ASSEMBLY_DIR)/assembly/total.fq $(ASSEMBLY_DIR)/assembly/single.fa
	mkdir -p $(ASSEMBLY_DIR)
	Trinity --seqType fq --JM 20G --single $(ASSEMBLY_DIR)/assembly/total.fq --CPU $(NUM_CPU) --output $(ASSEMBLY_DIR)/assembly 	

# Run gsAssembler on the basis of the SFF files in LIBRARY_SPEC
# note that the -ig -it and -icc options are trying to max out the algorithms effort
$(ASSEMBLY_DIR)/assembly/454Isotigs.fna: $(SOURCE_SFF_FILES) 
	newAssembly -cdna $(ASSEMBLY_DIR) 
	addRun -p $(ASSEMBLY_DIR) $(SOURCE_SFF_FILES)
	runProject -cpu $(NUM_CPU) -rip -acedir -m -ml $(MIN_LEN) -mi $(MIN_ID) -ig 10000 -it 10000 -icc 200 $(ASSEMBLY_DIR)

#gzip something
%.gz: %
	@ echo Compressing $< with gzip
	gzip $<

# Bowtie2 build
# there are actually a bunch of .bt2 files...
$(foreach i,1 2 3 4,%.$i.bt2): %
	$(BOWTIE2BUILD) $< $<

# Use bowtie2 to map an inputFile from LIBRARY_SPEC onto the cleaned-up OTUs
#$(SAM_DIR)/%.sam: $(FASTQ_DIR)/% $(foreach i,1 2 3 4,$(ASSEMBLED_OTU).seqclean.cd-hit.clstr.$i.bt2)
$(SAM_DIR)/%.sam: $(FASTQ_DIR)/% $(foreach i,1 2 3 4,$(ASSEMBLED_OTU).seqclean.cd-hit.$i.bt2)
	mkdir -p $(SAM_DIR)
#	$(BOWTIE2) --local -t -p $(NUM_CPU) $(ASSEMBLED_OTU).seqclean.cd-hit.clstr -U $< -S $@ --no-head
#	$(BOWTIE2) --local -t -p $(NUM_CPU) $(ASSEMBLED_OTU).seqclean.cd-hit -U $< -S $@ --no-head
	$(BOWTIE2) --local -t -p $(NUM_CPU) -x $(ASSEMBLED_OTU).seqclean.cd-hit -U $< -S $@ --no-head

# Create a CDB index for yanking
%.cidx: %
	$(CDBFASTA) $< -o $@

# Run SeqClean and clean up the litter
# see note above about how to increas SeqClean's cpu usage above 16
%.seqclean: %
	$(SEQCLEAN) $< -c $(NUM_CPU) -n 1 -v $(PRIMERS) -o $<.seqclean -r $<.cln
	rm -rf $(CWD)/cleaning_[1-9]* $(CWD)/outparts_cln.sort $(CWD)/err_seqcl_*.log $(CWD)/seqcl_*.log
# for some stoopid reason seqclean is dropping the .cidx file in the CWD
	rm $(CWD)/$(<F).cidx

# Winnow a FASTA file to remove the non-cpn60UT based on at least a 55% identity to anything in cpnDBNR
%.winnowed: % %.wateredBLAST
	$(mPUMA)/bin/winnow.pl -i $< -w $<.wateredBLAST -o $<.winnowed -c $(WINNOW_CUTOFF) 

# C2E is the portion of C3 that gets 150bps from each end of an OTU
%.c2e: %
	$(mPUMA)/bin/c2e.pl -i $< -o $<.c2e -l 150

# wateredBLAST is used to match each end against CPNDB_NR 
# the key idea in C3 is that if the ends don't seem to hit similar things that 
# its a putative chimera
# This also shows use of a simple FASTA batch processing trick to get some parallel execution going
%.wateredBLAST: %
# Call the watered BLAST on the ends of the OTUs in parallel based on the value of NUM_CPU
# N.B. this is actually untested on NUM_CPU == 1 but it should work
	$(mPUMA)/bin/fasta_batcher.pl -i $< -ni '-i' -o $<.wateredBLAST \
		-c $(NUM_CPU) -t $(NUM_CPU) -u 1 \
		-e '$(mPUMA)/bin/watered_blast.pl -d $(CPNDB_NR) -p blastn -h 0 -v 1'
	rm -f $(CWD)/fastaPiece*

# make the C3 report
%.c3report: % %.c2e.wateredBLAST
	$(mPUMA)/bin/c3.pl -i $< -w $<.c2e.wateredBLAST -h 0 -o $<.c3report

# make a list of the IDs that according to C3 are putatively chimeric
%.chimericIDs: %
	$(mPUMA)/bin/c3extractChimericIDs.pl -i $< -o $<.chimericIDs

# exclude the possible chimeras from the OTU list
%.possibleChimerasRemoved: % %.c3report.chimericIDs
	$(mPUMA)/bin/make_fasta_excluding_IDs.pl -i $< -e $<.c3report.chimericIDs -o $<.possibleChimerasRemoved 

# CDHit-EST @ 100% identity
# N.B. this is only for trivial condensation of OTUs. mPUMA strives to preserve as little as a single bp
$(OTU_PRE_CDHIT).cd-hit.clstr $(OTU_PRE_CDHIT).cd-hit: $(OTU_PRE_CDHIT)
	$(CDHITEST) -c 1 -G 1 -n 10 -M 0 -i $< -o $<.cd-hit

# need the CDHit at 100% identity for the amino acid space here
$(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit.clstr $(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit: $(OTU_PRE_CDHIT).cd-hit.oriented.aa
	$(CDHIT) -c 1 -G 1 -n 5 -M 0 -i $< -o $<.cd-hit

# need the FastTree for nucleotide
$(OTU_PRE_CDHIT).cd-hit.FastTree: $(OTU_PRE_CDHIT).cd-hit.aln
	$(FASTTREE) -nt -gtr $< > $@

# need the FastTree for the amino acid sequence
$(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit.FastTree: $(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit.aln
	$(FASTTREE) $< > $@

## There is a vestage in that mPUMA historically looks for a .fna file so just put down a symlink to support that
#%.cd-hit.clstr.fna: %.cd-hit 
#	ln -s $< $<.clstr.fna

# Produce a report from the CDHit-EST data 
# this can be used to see how OTUs should be collapsed on the basis of the clustering 100% identity in the DNA and
# AA sequence spaces
%.report: %
	$(mPUMA)/bin/parse_cd-hit-est.pl -c $< -o $<.report


# The major product of mPUMA is a table-of-contents (toc) file that shows how each read relates to the assembled
# OTUs
$(ASSEMBLY_DIR)/assembly/total.toc: $(ASSEMBLED_OTU) $(ASSEMBLY_DIR)/assembly/total.fna $(READ_TOC) \
		$(ASSEMBLED_OTU).seqclean $(OTU_PRE_CDHIT).cd-hit.clstr.report \
		$(ASSEMBLED_OTU).seqclean.c3report.chimericIDs 
	$(mPUMA)/bin/make_total_toc.pl -a $(ASSEMBLY_DIR)/assembly -total $(ASSEMBLY_DIR)/assembly/total.fna \
		-toc $(READ_TOC) -seqclean $(ASSEMBLED_OTU).seqclean \
		-clusterReport $(OTU_PRE_CDHIT).cd-hit.clstr.report \
		-chimericIDsFile $(ASSEMBLED_OTU).seqclean.c3report.chimericIDs \
		-useChimeraCalls $(USEC3)

# build the biom file
# N.B. that this will abide by WINNOW in terms of ensuring the reporting of only OTUs with sequences (i.e. winnowed set) to be in output
$(ASSEMBLY_DIR)/assembly/total.toc.biom: $(FINAL_OTU_SET) $(OTU_PRE_CDHIT).cd-hit.FastTree $(ASSEMBLY_DIR)/assembly/total.toc 
	$(mPUMA)/bin/mPUMA-to-biom-makefile.pl -a $(ASSEMBLY_METHOD) -m $(MAPPING_METHOD) -c $(USEC3) \
		-g $(GENENAME) -of $(FINAL_OTU_SET) -ff 1\
		-ot $(OTU_PRE_CDHIT).cd-hit.FastTree \
		-p $(ISPROTEIN) -o $@ -pp $(PRETTYPRINT) -v $(VERBOSE) -ls $(LIBRARY_SPEC) -toc $(ASSEMBLY_DIR)/assembly/total.toc 

# BLASTX of things against cpnDB-NR. This is used for determining orientation and then translating things
# formatdbPep($mPUMA::cpnDB_nr_pep);
%.blastx: % $(foreach i,pin psq phr, $(CPNDB_NR_PEP).$i)
	$(BLASTALL) -p blastx -i $< -o $<.blastx -d $(CPNDB_NR_PEP) -F F -v 100 -b 100 -a $(NUM_CPU)

# format a protein blast database
$(foreach i,pin psq phr,%.$i): %
	formatdb -i $< -p T

# format a nucleotide blast database
$(foreach i,nhr nin nsq,%.$i): %
	formatdb -i $< -p F

# Translate the proteins 
%.oriented.aa %.oriented: % %.blastx
	$(mPUMA)/bin/translate_stranded.pl -i $< -b $<.blastx -n $<.oriented -o $<.oriented.aa 

# t_coffee multiple sequence alignment 
%.aln: %
	$(TCOFFEE) $< -outfile $<.aln -mode quickaln -output fasta_aln

# clean up routines...
clean: clean_assembly

clean_assembly:
	rm -rf $(ASSEMBLY_DIR)

clean_condense:
	rm -f $(OTU_PRE_CDHIT).*

clean_abundance:
	rm -f $(READ_TOC) $(ALL_SAM_FILES) $(ASSEMBLY_DIR)/assembly/ace.toc

clean_biom:
	rm -f $(ASSEMBLY_DIR)/assembly/total.toc.biom

# stop just after the de novo assembly
assemble: $(ASSEMBLED_OTU)

# stop after condensing the OTUs
# the goal of OTU condensation is the SeqClean away primers, assess possible chimeras, 
# and then collapse all cases where a short OTU is 100% identical to a longer one
# also make the protein sequences
#condense: $(foreach i,1 2 3 4,$(OTU_PRE_CDHIT).cd-hit.clstr.$i.bt2) $(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit
condense: $(foreach i,1 2 3 4,$(OTU_PRE_CDHIT).cd-hit.$i.bt2) $(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit

# stop after the abundance data is completed
abundance: $(ASSEMBLY_DIR)/assembly/total.toc

# use tCoffee and FastTree for the phylogenetic workup (nucleotide and protein)
phylogeny: $(OTU_PRE_CDHIT).cd-hit.FastTree $(OTU_PRE_CDHIT).cd-hit.oriented.aa.cd-hit.FastTree

biom: $(ASSEMBLY_DIR)/assembly/total.toc.biom

PHONY: assemble condense abundance phylogeny biom all clean clean_assembly clean_condense clean_assembly

# future target for biom file?

