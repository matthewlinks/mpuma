#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use lib "$ENV{'mPUMA'}/lib";	# this allows us to find the module using the lone environment variable we require to be configured by the user
use mPUMA;

# Pipeline parameters: Initializing Assembly/Analysis Processes
my ($working_dir,$library_file,$name_of_assembly,$VERBOSE,$assembler);
my @options = (
	'w=s',	\$working_dir,
	'a=s',	\$assembler,
	'n=s',	\$name_of_assembly,
	'l=s',	\$library_file,	# This would be the 2-column text input which maps SFFs to libraries 
	'v=s',	\$VERBOSE
);
&GetOptions(@options);
usage() unless (
	defined $library_file
	and defined $name_of_assembly
);
$assembler = 'gsAssembler' unless defined $assembler;	# Default is gsAssembler / newbler

$working_dir = getcwd() unless defined $working_dir;
$mPUMA::VERBOSE = 1 if(defined($VERBOSE) eq 1);
mPUMA::checkWorkingDir($working_dir);

die "\$ENV{'mPUMA'} environment variable is undefined.\nPlease add the following command to the /home/username/.bashrc file\nexport mPUMA=/home/username/meta-cpn/\n\$ENV{'mPUMA'}" 
	unless defined $ENV{'mPUMA'};

# Parse the library file to figure out which NGS filess correspond to which libraries
my $library_info = mPUMA::parse_library_file( $library_file);	# parse_library_info must ensure no non-word chars in the library name
my @NGS_file_list;
foreach my $lib (keys %{$library_info}){
	push (@NGS_file_list,@{$library_info->{$lib}});
}

# make sure that the assembly routine is supported
if (defined $assembler){
	# do the assembly
	if ($assembler eq 'gsAssembler'){
		&assemble_using_gsAssembler($working_dir,$name_of_assembly,\@NGS_file_list) or die "Error calling assembly method: gsAssembler: $!";
	}elsif($assembler eq 'Trinity'){
		&assemble_using_Trinity($working_dir,$name_of_assembly,\@NGS_file_list) or die "Error calling assembly method: Trinity: $!";
	}else{
		die "Error assembly method $assembler is not supported";
	}
}else{
	die "Error there was no assembler defined";
}
exit 0;
sub usage {

die <<"USAGE";

	Usage: $0 -w working_dir -n name_for_assembly -l library_file -v VERBOSE -a gsAssembler

USAGE
}
sub assemble_using_gsAssembler{
	my $working_dir = shift;
	die "Error lost working directory" unless defined $working_dir;
	my $name_of_assembly = shift;
	die "Error lost assembly name" unless defined $name_of_assembly;
	my $NGS_file_list = shift;
	die "Error lost the list of NGS files to use" unless defined $NGS_file_list;


	# Initiates the project assembly process. Checks the outcome of the project assembly for any errors and troubleshoots any issues associated with the error.
	my $newblerDir = join('/',$working_dir,$name_of_assembly);	# what is the assembly going to be called
	my $assemblyDir = mPUMA::runAssembly( $newblerDir, @{$NGS_file_list});	# do the assembly and then return the path to the 'assembly' dir which newbler creates
}

sub assemble_using_Trinity{
	my $working_dir = shift;
	die "Error lost working directory" unless defined $working_dir;
	my $name_of_assembly = shift;
	die "Error lost assembly name" unless defined $name_of_assembly;
	my $NGS_file_list = shift;
	die "Error lost the list of NGS files to use" unless defined $NGS_file_list;

	my $trinityDir = join('/',$working_dir,$name_of_assembly);	# what is the assembly going to be called
	my $assemblyDir = mPUMA::runTrinity($trinityDir, @{$NGS_file_list});	# do the assembly and then return the path to the 'assembly' dir which is created

}
