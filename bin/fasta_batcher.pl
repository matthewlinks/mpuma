#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use IO::File;
use Parallel::Loops;
use Sys::CPU;
use Getopt::Long;
use File::Temp;

my ($fasta_file,$num_chunks,$num_threads,$output,$exe_cmd,$UNLINK,$VERBOSE,$name_input_param,$name_output_param);
my @options = (
        'i=s',  \$fasta_file,
	'ni=s',	\$name_input_param,
	'no=s',	\$name_output_param,
	'c=s',	\$num_chunks,
        't=s',	\$num_threads,
        'o=s',  \$output,
	'e=s',	\$exe_cmd,
	'u=s',	\$UNLINK,
	'v=s',	\$VERBOSE
);
&GetOptions(@options);

sub usage{
	die "Usage: $0 -i input_fasta -o output_to_create -e \'cmd with args\' -u 1 -c number_of_chunks_default_num_cpus -t number_of_threads_default_num_cpus";
}

usage() unless(
	defined $fasta_file
	and defined $exe_cmd
);

# Some defaults
$num_threads = Sys::CPU::cpu_count() unless defined $num_threads;
warn "Number of threads = $num_threads" if $VERBOSE;
$num_chunks = $num_threads unless defined $num_chunks;
warn "Number of chunks = $num_chunks" if $VERBOSE;
$UNLINK = 1 unless defined $UNLINK;
$name_input_param = '-i' unless defined $name_input_param;
$name_output_param = '-o' unless defined $name_output_param;

# check for obvious sheel no nos 
# N.B. this is only going to catch the most obvious errors
simple_sanity_check($exe_cmd);

# partition the input
warn "Splitting input into $num_chunks pieces" if $VERBOSE;
my $files = split_into_chunks($fasta_file,$num_chunks);

# create output filenames for the pieces 
# this needs to be in the scope of the master thread
my %input2output;
warn "Making tmpfiles for output" if $VERBOSE;
foreach my $input (@{$files}){
	my $fh = File::Temp->new(TEMPLATE => $input . '-output-XXXXX');
	my $output = $fh->filename;
	$input2output{$input} = $output;
}

# execute in parallel
my $parallel = Parallel::Loops->new($num_threads);
$parallel->share(\%input2output);

# run these in parallel
warn "Beginning parallel execution" if $VERBOSE;
$parallel->foreach($files,
	sub {
		my $input = $_;
		my $output = $input2output{$input};

		my $cmd = join(' ',$exe_cmd,
			$name_input_param,	$input,
			$name_output_param,	$output);

		warn "$$: Trying to run $cmd" if $VERBOSE;
		system($cmd) == 0 or die "Error calling $cmd: $!";
		warn "$$: has completed" if $VERBOSE;
		if ($VERBOSE){
			if (-e $output){
				warn "Child has tested and -e $output is true";
			}else{
				warn "Child has tested and -e $output is false";
			}
	
			if (-s $output){
				warn "Child has tested and -s $output is true";
			}else{
				warn "Child has tested and -s $output is false";
			}
		}
	}
);

# merge the outputs and return the results as requested 
# N.B. this does the unlink if $UNLINK is true
warn "Combining outputs" if $VERBOSE;
combine_outputs(\%input2output,$output);

undef $parallel;

exit 0;

sub simple_sanity_check {
	my $cmd = shift;
	die "Error lost command" unless defined $cmd;

	# check fot redictection and pipes
	if ($cmd =~ /\|/){
		die "Error $cmd contains a pipe";
	}elsif($cmd =~ /\`/){
		die "Error $cmd contains backticks";
	}elsif($cmd =~ /[\>\<]+/){
		die "Error $cmd contains redirection";
	}elsif($cmd =~ /-i\s+/){
		die "Error $cmd looks like it has an input parameter -i specified";
	}elsif($cmd =~ /-o\s+/){
		die "Error $cmd looks like it has an output parameter -o specified";
	}

	# check for likely input or output params -i and -o

}

sub combine_outputs{
	my $input2output = shift;
	die "Error lost structure with file names in it" unless defined $input2output;

	# write to a file if it was specified
	my $output = shift;

	warn "OUTPUT = [$output]" if $VERBOSE;

	if (defined $output) {
		open(STDOUT,">$output") or die "Error opening $output for writting: $!";
	}
	
	# get each of the sub results
	foreach my $input (keys %{$input2output}){
		my $out = $input2output->{$input};

		warn "Trying to read from output for $input from < $out" if $VERBOSE;
		open(OUT,"<$out") or die "Error opening $out for reading: $!";
		print STDOUT while(<OUT>);
		close OUT or die "Error closing $out: $!";
		# remove the output file we just read
		unlink($out) or die "Error unlinking $out: $!" if $UNLINK;
		# remove the input file 
		unlink($input) or die "Error unlinking $input: $!" if $UNLINK;
	}
	close STDOUT or die "Error closing STDOUT: $!";
}

sub unlink_files {
	my $input2output = shift;
	die "Error lost structure with file names in it" unless defined $input2output;
	foreach my $input (keys %{$input2output}){
		unlink($input2output->{$input}) or die "Error unlinking ",$input2output->{$input},": $!";
		unlink($input) or die "Error unlinking ",$input,": $!";
	}
}

sub split_into_chunks {
	my $file = shift;
	die "Error lost file" unless defined $file;
	my $n_chunks = shift;
	die "Error lost the number of chunks" unless defined $n_chunks;

	# get a bunch of temp files
	my $files = get_files($n_chunks);
	die "Error lost output streams" unless defined $files;
	my $num_files = @{$files};
	die "Error wrong number of files $num_files" unless ($num_files == $n_chunks);

	# open each file as a Bio::SeqIO
	my @fastaIOs;
	foreach my $f (@{$files}){
		warn "Writting to $f" if $VERBOSE;
		my $seqio = Bio::SeqIO->new(-file => ">$f", -format => 'fasta');
		push(@fastaIOs,$seqio);
	}

	# split the input file
	my $inio = Bio::SeqIO->new(-file => "<$file", -format => 'fasta');
	my $i = 0;
	while(my $seq = $inio->next_seq()){
		$fastaIOs[$i]->write_seq($seq);
		$i++;
		$i = 0 unless ($i < @fastaIOs);
	}
	foreach my $seqio (@fastaIOs){
		undef $seqio;
	}

	return $files;
}

sub get_files {
	my $num = shift;
	die "Error lost num" unless defined $num;

	my @files;
	for(my $i = 1; $i <= $num; $i++){
		my $fh = File::Temp->new(TEMPLATE => 'fastaPieceXXXXX');
		push(@files,$fh->filename);
	}
	return \@files;
}
