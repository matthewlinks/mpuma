#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use lib "$ENV{'mPUMA'}/lib";    # this allows us to find the module using the lone environment variable we require to be configured by the user
use mPUMA;
my ($assemblyDir,$totalReadFile,$aceToc,$seqcleanFile,$clusterReport,$chimericIDsFile,$useChimeraCalls);
my @options = (
	'a=s',			\$assemblyDir,
	'total=s',		\$totalReadFile,
	'toc=s',		\$aceToc,
	'seqclean=s',		\$seqcleanFile,
	'clusterReport=s',	\$clusterReport,
	'chimericIDsFile=s',	\$chimericIDsFile,
	'useChimeraCalls=s',	\$useChimeraCalls	
);
&GetOptions(@options);

die "Usage: $0 -a assemblyDir -total totalReadFile -toc sam.toc -seqclean something.seqclean -clusterReport something.cd-hit.clstr.report -chimericIDsFile something.chimericIDs =useChimeraCalls 0"
unless defined $assemblyDir and defined $totalReadFile and defined $aceToc and defined $seqcleanFile and defined $clusterReport and defined $chimericIDsFile and defined $useChimeraCalls;

mPUMA::make_total_toc($assemblyDir,$totalReadFile,$aceToc,$seqcleanFile,$clusterReport,$chimericIDsFile,$useChimeraCalls);
