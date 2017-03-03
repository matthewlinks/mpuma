#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($assemblyDir,$database);
my @options = (
	'a=s', \$assemblyDir,
	'd=s',	\$database
);
&GetOptions(@options);
die "Usage: $0 -a Newbler_assembly_directory -d database" unless (
	defined $assemblyDir and defined $database);

make_total_msf($assemblyDir);
sub make_total_msf {
	my $dir = shift;
	die "Error lost directory unless defined $dir: $!";
	my $total_msf = join('/',$dir,'total.fna');
	die "Error $total_msf already exists" if (-e $total_msf);
	my $sff_dir = join('/',$dir,'sff');

	open(SFFDIR,join


# run the assembly 
time runProject -cpu 10 -rip -acedir -m -ml 137 -mi 90 -vt /scratch/linksm/MVP/cpn60-primers-enumerated.msf

# waterBLAST the isotigs
$METACPN/bin/watered_blast.pl -v 5 -p blastn -i 454Isotigs.fna -d /aped/blast_dbs/cpndb_nr_20100504 > 454Isotigs.fna.waterBLAST.cpndb_nr_20100504

# create the total read set
# something like for i in sff/*.sff; do sffinfo -s $i >>total.fna; done
for i in `cat SFF.list`; do echo $i; sffinfo -s /aped/sffs/$i >> total.fna; done

# waterBLAST the individual reads
# this is required to check the Specificity and Sensitivity of assignments
$METACPN/bin/watered_blast.pl -v 5 -p blastn -i total.fna -d /aped/blast_dbs/cpndb_nr_20100504 > total.fna.waterBLAST.cpndb_nr_20100504
# Atul's data
#real    224m18.478s
#user    205m40.636s
#sys     67m46.555s
# Bonnie 
#real    191m41.398s
#user    173m48.888s
#sys     53m26.837s

# calculate the Specificity and Sensitivity metrics 
$METACPN/bin/calculate_Sp_Sn_metrics.pl -a /path/to/assembly/ -i /path/to/total.fna -v 5 -d /aped/blast_dbs/cpndb_nr_20100504 -verbose 1 -s 1
#time $METACPN/bin/calculate_Sp_Sn_metrics.pl -a /home/pipeline/Atul-Bonnie-pyro/assemblies/Atul-trout-July2010/assembly/ -i /home/pipeline/Atul-Bonnie-pyro/assemblies/Atul-trout-July2010/total.fna -v 5 -d /aped/blast_dbs/cpndb_nr_20100504 -verbose 1 -s 1 > Specificity_Sensitivity

# create directories 1 / project
# each directory is a library we want to analyse / compare
...

# create a file called files in each project directory which lists the SFF files which make up that library
...

# create a reads.ids file in each project directory which lists the read specific to that library
for i in *; do echo $i; for j in `cat $i/files`; do sffinfo -a /aped/sffs/$j >> $i/reads.ids; done; done

# create the diversity files for each project

time for i in *; do echo $i; perl $METACPN/bin/generate_diversity_subset.pl -a /path/to/assembly/ -aw /path/to/assembly/454Isotigs.fna.wateredBLASTN -i /path/to/total.fna -iw /path/to/total.fna.wateredBLASTN -s $i/reads.ids -scale 20000 > $i/diversity.scaled20k.txt; done
#time for i in BC-*; do echo $i; perl $METACPN/bin/generate_diversity_subset.pl -a /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010/assembly/ -aw /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010/assembly/454Isotigs.fna.wateredBLASTN -i /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010/total.fna -iw /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010/total.fna.wateredBLASTN -s $i/reads.ids -scale 4491 > $i/diversity.scaled4491.txt; done
#time for i in AD-*; do echo $i; perl $METACPN/bin/generate_diversity_subset.pl -a /home/pipeline/Atul-Bonnie-pyro/assemblies/Atul-trout-July2010/assembly/ -aw /home/pipeline/Atul-Bonnie-pyro/assemblies/Atul-trout-July2010/assembly/454Isotigs.fna.wateredBLASTN -i /home/pipeline/Atul-Bonnie-pyro/assemblies/Atul-trout-July2010/total.fna -iw /home/pipeline/Atul-Bonnie-pyro/assemblies/Atul-trout-July2010/total.fna.wateredBLASTN -s $i/reads.ids -scale 1436 > $i/diversity.scaled1436.txt; done

# Updates for the -rip assembly
#for i in `cat ../BC-normal` ; do echo $i; perl $METACPN/bin/generate_diversity_subset.pl -a /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/assembly/ -aw /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/assembly/454Isotigs.fna.wateredBLASTN  -s BC-normal/$i/reads.ids -scale 6021 -i /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/total.fna -iw /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/total.fna.wateredBLASTN > BC-normal/$i/diversity.scaled.6021.txt; done

#for i in `cat ../BC-enriched` ; do echo $i; perl $METACPN/bin/generate_diversity_subset.pl -a /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/assembly/ -aw /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/assembly/454Isotigs.fna.wateredBLASTN  -s BC-enriched/$i/reads.ids -scale 1792 -i /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/total.fna -iw /home/pipeline/Atul-Bonnie-pyro/assemblies/Bonnie-dog-July2010-rip/total.fna.wateredBLASTN > BC-enriched/$i/diversity.scaled1792.txt; done

# create the genespring files / library
#for i in `cat ../BC-normal` ; do echo $i; grep -e "^isotig" BC-normal/$i/diversity.scaled.6021.txt | sort > ./for-genespring-BC-normal/$i.txt; done

# create an isotig-list file
#cat *.txt | grep -e "isotig" | cut -f 1  | sort | uniq > isotig-list

# create a tech file which needs the header deleted...or a copy created
#echo -e "OTU\tPercent\tScaled6021\tActualCount\tLabel" > technology.tsv; for i in `cat isotig-list`; do echo -e "$i\t0\t0\t0\tsomething"; done >> technology.tsv

# convert the technology around so that it contains a row for each OTU
#/home/pipeline/meta-cpn/bin/convert_technology.pl -t technology.tsv -o ../ *.txt

C'est ca
