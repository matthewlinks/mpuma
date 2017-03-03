#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use lib "$ENV{'mPUMA'}/lib";    # this allows us to find the module using the lone environment variable we require to be configured by the user
use mPUMA;

my ($assemblyDir,$totalReadFile,$samToc,$seqcleanFile,$clusterReport,$chimericIDsFile,$useChimeraInfo,$VERBOSE);

my @options = (
	'a=s',	\$assemblyDir,
	't=s',	\$totalReadFile,
	's=s',	\$samToc,
	'sc=s',	\$seqcleanFile,
	'cr=s',	\$clusterReport,
	'ci=s',	\$chimericIDsFile,
	'uc=s',	\$useChimeraInfo,
	'v=s',	\$VERBOSE);
&GetOptions(@options);

die "Usage: $0 -a ./ -t total.fna -s sam.toc -sc Trinity.fasta.seqclean -cr Trinity.fasta.seqclean.cd-hit.clstr.report -ci Trinity.fasta.seqclean.c3report.chimericIDs -uc 0"
unless (
	defined $assemblyDir
	and defined $totalReadFile
	and defined $samToc
	and defined $seqcleanFile
	and defined $clusterReport
	and defined $chimericIDsFile
	and defined $useChimeraInfo
);

$mPUMA::VERBOSE = $VERBOSE;

mPUMA::make_total_toc($assemblyDir,$totalReadFile,$samToc,$seqcleanFile,$clusterReport,$chimericIDsFile,$useChimeraInfo);


exit 0;
