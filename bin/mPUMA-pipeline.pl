#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use lib "$ENV{'mPUMA'}/lib";	# this allows us to find the module using the lone environment variable we require to be configured by the user
use mPUMA;
my ($blast_db,$working_dir,$library_file,$name_of_assembly,$VERBOSE,$CHIMERACHECK,$ASSEMBLER,$MAPPING);
my @options = (
	'a=s',	\$ASSEMBLER,
	'mapping=s',	\$MAPPING,
	'w=s',	\$working_dir,
	'n=s',	\$name_of_assembly,
	'l=s',	\$library_file,	# This would be the 2-column tet input which maps SFFs to libraries 
	'v=s',	\$VERBOSE,
	'c=s',	\$CHIMERACHECK,
);
&GetOptions(@options);
usage() unless (
	defined $library_file
	and defined $working_dir
);

$VERBOSE	= 0 unless defined $VERBOSE;
$CHIMERACHECK	= 1 unless defined $CHIMERACHECK; 
$ASSEMBLER	= 'gsAssembler' unless defined $ASSEMBLER;
if (!defined $MAPPING){
	if ($ASSEMBLER eq 'gsAssembler'){
		$MAPPING = 'bowtie2';
	}elsif ($ASSEMBLER eq 'Trinity'){
		$MAPPING = 'bowtie2' 
	}else{
		die "Error in unsupported assembler mode";
	}
}
mPUMA::checkWorkingDir($working_dir);

sub usage {

die <<"USAGE";

	Usage: $0 -w working_dir -l library_file

	This runs the mPUMA assembly process and the post assembly analysis with the default settings. 

	   OPTIONS:
	   -a Assembler		Default: gsAssembler optional: Trinity
	   -mapping method	Default: bowtie2 optional but not recommended: gsAssembler
	   -c boolean		Default: 1 - set to 0 to ignore C3 chimera report.
	   -n AssemblyName	Default: Timestamp
	   -v VERBOSE		Default: 0

USAGE
}

die "\$ENV{'mPUMA'} environment variable is undefined.\nPlease add the following command to the /home/username/.bashrc file\nexport mPUMA=/home/username/meta-cpn/\n\$ENV{'mPUMA'}" 
	unless defined $ENV{'mPUMA'};

$name_of_assembly = generate_assembly_name() unless (defined $name_of_assembly);

# this is basically a 2-step process. 
# Step 1 - generate an assembly
# Step 2 - convert the assembled position of reads into abundance information and output data in formats for other downstream tools

# Step 1 - call the assembly script
system("$ENV{mPUMA}/bin/mPUMA-pipeline-assembly.pl",
	'-v',	$VERBOSE,
	'-w',	$working_dir,
	'-n',	$name_of_assembly,
	'-a',	$ASSEMBLER,
	'-l',	$library_file) == 0 or die "Error calling mPUMA-pipeline-assembly.pl: $!";

# Step 2 - call the post-assembly analysis script
system("$ENV{mPUMA}/bin/mPUMA-pipeline-analysis.pl",
	'-v',		$VERBOSE,
	'-c',		$CHIMERACHECK,
	'-w',		$working_dir,
	'-n',		$name_of_assembly,
	'-a',		$ASSEMBLER,
	'-mapping',	$MAPPING,
	'-l',		$library_file) == 0 or die "Error calling mPUMA-pipeline-analysis.pl: $!";

# Give the user a nudge to try re-running mPUMA-pipeline-analysis.pl if they want to modify the outputs
print "\n\n\n==========DONE==========\nmPUMA has finished processing $name_of_assembly in $working_dir. If you want to manage additional outputs you may want to rerun \$mPUMA/bin/mPUMA-pipeline-analysis.pl modifying its parameters accordingly.\n";

exit 0;

sub generate_assembly_name {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $string = $dayOfMonth . $months[$month] . $year;
	return $string;
}
