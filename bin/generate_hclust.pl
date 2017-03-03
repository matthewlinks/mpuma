#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

# Program for the hierarchial clustering of OTUs between multiple libraries.

# Command line example for supplying program with the proper parameter options.
# perl $mPUMA/bin/generate_hclust.pl -i ~/tomato_analysis_20June2012_nway_comparisons/nucleotide_comparisons/pairwise_shared_comparisons/S1-TH-8_S2-TH-9G/PARSED_R_INPUT_FILES -o ~/HCLUST_FILES

# Declaration of parameter options to be initialized in @options argument strings.
my ($input_dir, $output_dir, $R_suffix, $hclust_suffix, $distFunction, $hclustFunction, $logbase, $scale, $dendrogram, $trace);
my @options = (
	'i=s',	\$input_dir,
	'o=s',	\$output_dir,
	's=s',	\$R_suffix,
	't=s',	\$hclust_suffix,
	'd=s',	\$distFunction,
	'h=s',	\$hclustFunction,
	'l=s',	\$logbase,
	'c=s',	\$scale,
	'g=s',	\$dendrogram,
	'r=s',	\$trace
);
&GetOptions(@options);

usage() unless (
	defined $input_dir
	and defined $output_dir
);
# Parameter options initialized if not defined in @options argument strings.

# R input file suffix. Can be anything (R.input.txt, etc.)
$R_suffix = 'R.input.txt' unless defined $R_suffix;

# HCLUST input file suffix. Can be anything (hclust.txt, etc.)
$hclust_suffix = 'hclust.txt' unless defined $hclust_suffix;

# The distance measure to be used on the data. This should be (an unambiguous abbreviation of or full word) one of either "euclidean", "maximum", "manhattan" or "binary". The default is "euclidean".
$distFunction = 'euclidean' unless defined $distFunction;

# The agglomeration method to be used in hierarchial clustering. This should be (an unambiguous abbreviation of or full word) one of either "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". The default is "average".
$hclustFunction = 'average' unless defined $hclustFunction;

# The logbase to perform a log transform on the datamatrix.
$logbase = 10 unless defined $logbase;

# Scale the values on the row and/or column of the infile. Can be specified as "row", "column", or "none". DEFAULT is "none".
$scale = "none" unless defined $scale;

# Draw dendrograms on either the row and/or column of the infile. Can be specified as "row", "column", "both", or "none". DEFAULT is "both".
$dendrogram = "both" unless defined $dendrogram;

# Trace the row and/or column of the infile. Can be specified as "row", "column", "both", or "none". DEFAULT is "none".
$trace = "none" unless defined $trace;

# Usage details on input parameters for the hierarchial clustering of shared OTUs between multiple libraries.
sub usage {

die <<"USAGE";

Usage: $0 -i input_dir -o output_dir -d distFunction -h hclustFunction -s R_suffix -t hclust_suffix -c scale -g dendrogram -r trace

OPTIONS:

-i input_dir - The input_dir containing files in the following format. Where the header consists of the libraries to compare and each subsequent line indicates the isotig and the scaled values found within the libraries specified in the header line.
e.g.
S1-TH-8-1	S1-TH-8-2	S1-TH-8-3	S1-TH-8-4	S1-TH-8-5	S1-TVD-1-1	S1-TVD-1-2	S1-TVD-1-3	S1-TVD-1-4	S1-TVD-1-5
isotig00001	16.00	45.00	121.00	125.00	40.00	43.00	229.00	53.00	65.00	96.00

-o output_dir - The directory to output the hclust associated images and files.

-s R_suffix - R input file suffix. Can be anything (R.input.txt, etc.)

-t hclust_suffix - Hclust input file suffix. Can be anything (hclust.txt, etc.)

-d distFunction - The distance measure to be used on the data. This should be (an unambiguous abbreviation of or full word) one of either "euclidean", "maximum", "manhattan" or "binary". The default is "euclidean".

-h hclustFunction - The agglomeration method to be used in hierarchial clustering. This should be (an unambiguous abbreviation of or full word) one of either "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". The default is "average".

-c scale - Scale the values on the row and/or column of the infile. Can be specified as "row", "column", or "none". DEFAULT is "none".

-d dendrogram - Draw dendrograms on either the row and/or column of the infile. Can be specified as "row", "column", "both", or "none". DEFAULT is "both".

-t trace - Trace the row and/or column of the infile. Can be specified as "row", "column", "both", or "none". DEFAULT is "none".

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the R input file directory for any R input files with the specified suffix.
my $R_input_files = find_input_files($input_dir,$R_suffix);

# Foreach input file in the R input file list perform hierarchial clustering on shared OTUs between multiple libraries.
foreach my $filename (@{$R_input_files}){
	warn "$filename\n";
	# Get the R input file path.
	my $R_infile = join('/', $input_dir, $filename);
	
	# Extract library names from R.input filename for the hclust dendrogram plot title and the outfile name.
	my ($title, $outfile);
	if($filename =~ m/([\w\_\-]+)\.$R_suffix/){
		$outfile = $1;
		# Remove shared tag
		$outfile =~ s/.shared//g;
		my @split_filename = split(/_/, $outfile);
		$title = join(" vs. ", @split_filename);
	}
	
	# The HCLUST libraries output filename.
	my $hclust_libraries_outfile = join('/', $output_dir, $outfile . "-libraries");

	# The HCLUST OTUs output filename.
	my $hclust_OTUs_outfile = join('/', $output_dir, $outfile . "-OTUs");

# Pipe the following R script that computes hierarchial clustering method on the R input file to the R program.
open (R_SCRIPT_PIPE, "|R --vanilla --slave");
print R_SCRIPT_PIPE <<EOF;
## Loading Global library packages.

## Cluster and Tree Conversion Package for writting Newick Tree files.
library(ctc)

## Read the data file into a data frame.
dataframe <- read.delim("$R_infile")

datacolumns <- names(dataframe)
datamatrix <- as.matrix(dataframe[datacolumns])
datamatrix[datamatrix==0] <- 1
datamatrix <- log(datamatrix, $logbase)

if (!file.exists("$hclust_libraries_outfile-dendrogram.png")){
	## Libraries are the direct matrix
	libraries <- (datamatrix)

	## Hierarchical Clustering
	## Generate the libraries HCLUST dendrogram.
	hclustLibraries <- hclust(dist(libraries, method = "$distFunction"), method="$hclustFunction")
	png('$hclust_libraries_outfile-hclust-dendrogram.png\');
	plot(hclustLibraries, hang = -1, main="Cluster Dendrogram ( $title )")
	
	## Generate Libraries table from HCLUST results.
	hclustLibrariesDendroLabels <- data.frame(labels=hclustLibraries\$labels[hclustLibraries\$order])
	write.table(hclustLibrariesDendroLabels, file = "$hclust_libraries_outfile.hclust.txt", sep = "\t")

	## Generate Libraries newick tree from HCLUST results.
	write(hc2Newick(hclustLibraries),file='$hclust_libraries_outfile-hclust.newick')

	## Turn off R graphical device
	## If successful, prints on STDOUT
	## null device
	##          1
	dev.off()

} else {
	print("$outfile-libraries-hclust-dendrogram.png has been previously generated so skipping...")
}
# if (!file.exists("$hclust_OTUs_outfile-hclust-dendrogram.png")){
if (!file.exists("$hclust_OTUs_outfile-hclust.newick")){
	## OTUs are the transposed matrix
	OTUs <- t(datamatrix)

	## Hierarchical Clustering
	## Generate the OTUs HCLUST dendrogram.
	hclustOTUs <- hclust(dist(OTUs, method = "$distFunction"), method="$hclustFunction")
	png('$hclust_OTUs_outfile-hclust-dendrogram.png\');
	plot(hclustOTUs, hang = -1, main="Cluster Dendrogram ( $title )")

	## Generate OTUs table from HCLUST results.
	hclustOTUsDendroLabels <- data.frame(labels=rev(hclustOTUs\$labels[hclustOTUs\$order]))
	write.table(hclustOTUsDendroLabels, file = "$hclust_OTUs_outfile.hclust.txt", sep = "\t")

	## Generate OTUs newick tree from HCLUST results.
	write(hc2Newick(hclustOTUs),file='$hclust_OTUs_outfile-hclust.newick')

	## Turn off R graphical device
	## If successful, prints on STDOUT
	## null device
	##          1
	dev.off()

} else {
	print("$outfile-OTUs-hclust-dendrogram.png has been previously generated so skipping...")
}

## Exit out of R program
q()

EOF
close(R_SCRIPT_PIPE);


# The heatmap output image filename.
my $heatmap_outfile = join('/', $output_dir, $outfile . "-heatmap.pdf");

# Pipe the following R script to the R program that generates a heatmap using the R input file.
# open (R_SCRIPT_PIPE, "|R --vanilla --slave");
open (R_SCRIPT_PIPE, "|R --vanilla --slave");
print R_SCRIPT_PIPE <<EOF;
if (!file.exists("$heatmap_outfile")){
	## R Color Brewer Package
	library(RColorBrewer)
	## Pretty Heatmap Package
	library(pheatmap)
	
	## Read the data file into a data frame.
	datamatrix <- as.matrix(read.table("$R_infile", header = TRUE,sep="\t"))
	# edit for re-orienting the input
	datamatrix <- t(datamatrix)

	datamatrix[datamatrix==0] <- 1
	## Generate Heatmap
	pheatmap(log(datamatrix, $logbase),
 	filename="$heatmap_outfile",
	border_color = NA,
	cellwidth = 60,
	cellheight = 12,
	color=colorRampPalette(c("blue","yellow","red"))(256),
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	clustering_distance_cols = "$distFunction",
	clustering_distance_rows = "$distFunction",
	clustering_method = "$hclustFunction",
	bg = "transparent",
	treeheight_row = 0,
	treeheight_col = 500,
	legend = TRUE,
	)

	## Turn off R graphical device
	## If successful, prints on STDOUT
	dev.off()
	
	## Exit out of R program
	q()
} else {
	print("$outfile-heatmap.png has been previously generated so skipping...")
}
EOF
close(R_SCRIPT_PIPE);
	
}

# Parse the hclust input file directory for any hclust input files with the specified suffix.
my $hclust_input_files = find_input_files($output_dir, $hclust_suffix);

# Foreach input file in the hclust input file list parse the hclust isotig or sample lists for a cleaner file.
foreach my $filename (@{$hclust_input_files}){

	# Get the hclust input file path.
	my $hclust_infile = join('/', $output_dir, $filename);

	my $outfile;
	if($filename =~ m/([\w\_\-]+)\.$hclust_suffix/){
		$outfile = $1;
	}
	my $hclust_outfile = join('/', $output_dir, $outfile . "-hclust-dendrogram.txt");
	open INFILE, "<$hclust_infile" or die "Error opening $hclust_infile for reading: $!";
	open OUTFILE, ">$hclust_outfile" or die "Error opening $hclust_outfile for writting: $!";
	while(<INFILE>){
		chomp $_;
		if ($_ =~ m/^"\d+"\t"(isotig\d+|[\w\_\-\.]+)"$/){
			my $clustered_item = $1;
			$clustered_item =~ s/\./-/g;
			print OUTFILE "$clustered_item\n";
		}
		
	}
	close INFILE or die "Error closing $hclust_infile: $!";
	close OUTFILE or die "Error closing $hclust_outfile: $!";
	if(-s $hclust_outfile){
		unlink $hclust_infile;
	}
}



###
exit 0;

=head2 B<\@input_files> = B<find_input_files(B<$input_dir>, B<$suffix>)> - Find all input files with the directory and with the suffix specified and returns a reference to the list of input files

I<Input Parameters>:

B<$input_dir>  - The input file directory to parse the input files.

B<$suffix>  - Input file suffix. Can be anything (R.input.txt, etc.)

I<Output Parameters>:

B<\@input_files> - Reference to a list of input files.

E<10>

E<10>

=cut
sub find_input_files {

	# The directory to parse for input files.
	my $input_dir = shift;
	die "Error lost input file directory" unless defined $input_dir;

	# Input file suffix. Can be anything (R.input.txt, etc.)
	my $suffix = shift;
	die "Error lost input file suffix" unless defined $suffix;

	if(-d $input_dir){ # Check if $input_dir is a directory.
		opendir(DIR,$input_dir) or die "Error opening $input_dir: $!";
		my @input_files;
		while(my $input_file = readdir(DIR)){
			warn "Checking $input_file for $suffix";
			if ($input_file =~ /\.$suffix$/){
				push(@input_files, $input_file);
			}
		}
		closedir(DIR) or die "Error closing $input_dir: $!";

		# Reference to a list of input files.
		return \@input_files;
	}
	else{
		die "Error $input_dir does not exist!\n";
	}
}
