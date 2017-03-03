#!/usr/bin/perl -w
use strict;
use lib "$ENV{mPUMA}/lib";
use mPUMA;
use Getopt::Long;

my ($subDir);
my @options = (
	's=s',	\$subDir
);
&GetOptions(@options);
$subDir = './' unless defined $subDir;


my $sffFile = join('/',$subDir,'reads.sff');
die "Error $sffFile does not exist." unless (-e $sffFile);

my $timeStamp = generate_assembly_name();
my $dir = join('/',$subDir,$timeStamp);

my $assemblyDir = mPUMA::newAssembly($dir);
mPUMA::addRun($dir,$sffFile);
mPUMA::runProject($assemblyDir);

print "Check $dir\n";

exit 0;


sub generate_assembly_name {
        my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
        my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
        my $year = 1900 + $yearOffset;
        my $string = $dayOfMonth . $months[$month] . $year;
        return $string;
}
