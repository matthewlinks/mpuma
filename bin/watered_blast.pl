#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::Perl;
use Bio::SeqIO;
use Bio::SearchIO;
use IPC::Open2;
use IO::File;
use POSIX qw(tmpnam);

$| = 1;

my ($input_file,$blast_db,$blast_prog,$num_descs,$output_file,$print_header);
my @options = (
	'i=s',	\$input_file,
	'd=s',	\$blast_db,
	'p=s',	\$blast_prog,
	'v=s',	\$num_descs,
	'o=s',	\$output_file,
	'h=s',	\$print_header,
);
&GetOptions(@options);
$num_descs = 5 unless defined $num_descs;
usage() unless (
	defined $input_file
	and defined $blast_db
	and defined $blast_prog
);
$print_header = 0 unless defined $print_header;
die "Error BLAST program must be either blastp or blastn. This is because a Smith - Waterman alignment is performed on the top hit(s) and thus the sequences must be of the same type" 
	unless (($blast_prog eq 'blastp') or ($blast_prog eq 'blastn'));

# redirect stdout if there is a file to write to
open(STDOUT,">$output_file") or die "Error opening $output_file for writting: $!" if (defined $output_file);

# Ensure BLAST formated
ensure_blast_formated($blast_db,$blast_prog);

# N.B. cdbindexing is ensured in the routine which calles cdbyank so don't worry about it up here

# get a blastIO to the BLAST search
# N.B. remeber to cleanup with a wait later
my ($blastIO,$pid) = get_blast_io($input_file,$blast_db,$blast_prog,$num_descs);

print join("\t",'Query','Description','% ID of match (not to self)','Length','Closest hit header','Strandedness'),"\n" if ($print_header);;

# Walk through the Blast report printing the results on the limit of the sequence identity
while (my $result = $blastIO->next_result() ) {
	my $hit = $result->next_hit();
	if (!defined $hit){
#		warn "Error not a single hit for: ",$result->query_name();
		print join("\t",$result->query_name(),$result->query_description(),0,'N/A','N/A'),"\n";
		next;
	}

	my $best_significance = $hit->significance();

	# Make sure $hit points to the first different sequence name
	# this is checked for the case where we were given the same sequence for the FASTA and the BLAST database
	my $q_name = $result->query_name();
	my $h_name = $hit->name();
	until($h_name ne $q_name){
		$hit = $result->next_hit();
		next unless defined $hit;
		$h_name = $hit->name();
	}

	# Get the Bio::PrimarySeqI for the query and the hit
	my $q_seq = cdbyank_seq($result->query_name(),$input_file);
	my $h_seq = cdbyank_seq($hit->name(),$blast_db);

	# handle strandedness
	my $qstrand = $hit->strand('query');
	my $sstrand = $hit->strand('hit');
	my $strand = '+';
	if ($qstrand != $sstrand){
		$strand = '-';
	}


	# figure out what the % identity is to the next best hit
	my $limit_percent_identity = find_percent_id($q_seq,$h_seq,$strand);
	die "Error getting percent ID limit for ",$result->query_name()," <=> ",$hit->name() unless defined $limit_percent_identity;

	# print the results
	print join("\t",$q_seq->display_id(),$q_seq->desc(),$limit_percent_identity,$q_seq->length(),$h_seq->display_id() . ' ' . $h_seq->desc(),$strand),"\n";

	# add a line for each hit at the same significance
	while($hit = $result->next_hit()){
		last unless ($hit->significance() == $best_significance);
#		warn "Multiples for $q_name";
		$h_name = $hit->name();
		until($h_name ne $q_name){
			$hit = $result->next_hit();
			next unless defined $hit;
			$h_name = $hit->name();
		}
	
		# Get the Bio::PrimarySeqI for the query and the hit
		my $q_seq = cdbyank_seq($result->query_name(),$input_file);
		my $h_seq = cdbyank_seq($hit->name(),$blast_db);
		my $qstrand = $hit->strand('query');
		my $sstrand = $hit->strand('hit');
		my $strand = '+';
		if ($qstrand != $sstrand){
#			warn "On $q_name <=> $h_name got strands $qstrand $sstrand";
			$strand = '-';
		}
	
		# figure out what the % identity is to the next best hit
		my $limit_percent_identity = find_percent_id($q_seq,$h_seq,$strand);
		die "Error getting percent ID limit for ",$result->query_name()," <=> ",$hit->name() unless defined $limit_percent_identity;
	
		# print the results
		print join("\t",$q_seq->display_id(),$q_seq->desc(),$limit_percent_identity,$q_seq->length(),$h_seq->display_id() . ' ' . $h_seq->desc(),$strand),"\n";
	}		
}
cleanup();

exit 0;

sub find_percent_id{
	my $q_seq = shift;
	die "Error lost query sequence" unless defined $q_seq;
	my $h_seq = shift;
	die "Error lost hit sequence" unless defined $h_seq;
	my $qstrand = shift;
	die "Error lost strandedness of the query" unless	


	# Write the files to temporary files
	my $q_file = write_seq_to_tmp_file($q_seq,$qstrand);
	my $h_file = write_seq_to_tmp_file($h_seq,'+');

	# call water 
	my $pid = open2(\*WATER_OUT,\*WATER_IN,'water',
		'-asequence',	$q_file,
		'-bsequence',	$h_file,
		'-gapopen',	'10.0',
		'-gapextend',	'0.5',
		'-stdout',	'-auto') or die "Error calling water: $?";
	close WATER_IN or die "Error closing STDIN to water process: $!";
#	warn "Water call water -asequence $q_file -bsequence $h_file -gapopen 10.0 -gapextend 0.5 -stdout -auto";

	my $id;
	while(<WATER_OUT>){
		chomp $_;
		if ($_ =~ /^#\s+Identity\:\s+\d+\/\d+\s+\(([\d\.]+)\%\)/){
			$id = $1;
			last;
		}
	}
	close WATER_OUT or die "Error closing STDOUT from water process: $!";
	wait;

	# cleanup the tmp files
	unlink($q_file) or die "Error unlinking $q_file: $!";
	unlink($h_file) or die "Error unlinking $h_file: $!";

	return $id;
}

sub write_seq_to_tmp_file {
	my $seq = shift;
	die "Error lost sequence" unless defined $seq;
	my $strand = shift;
	die "Error lost strand" unless defined $strand;

	# reverse complement as necessary
	if ($strand eq '-'){
		$seq->seq(reverse_complement_as_string($seq->seq()));
	}elsif($strand ne '+'){
		die "Unknown strand $strand";
	}

	my $name = get_tmp_file();
	my $io = Bio::SeqIO->new(-file => ">$name", -format => 'fasta');
	$io->write_seq($seq);
	return $name;
}
	
sub get_tmp_file {
	my $name;	
	do { $name = tmpnam() }
    	until my $fh = IO::File->new($name, O_RDWR|O_CREAT|O_EXCL);
	# N.B. The calling process is responsible for unlinking this file!!!
	return $name;
}

sub cdbyank_seq{
	my $id = shift;
	die "Error lost id" unless defined $id;
	my $fasta_file = shift;
	die "Error lost fasta file" unless defined $fasta_file;
	my $cdbfasta_file = $fasta_file . '.cidx';
	if (!(-e $cdbfasta_file)){
		die "Error fasta file does not exist" unless (-e $fasta_file);
		# call cdbfasta on the input file if the indexed file does not exist already
		system('cdbfasta',$fasta_file,'-o', $cdbfasta_file) == 0 or die "Error calling cdbfasta $cdbfasta_file: $?";
	}
        my $pid = open2(\*CHLD_OUT, \*CHLD_IN,'cdbyank',$cdbfasta_file) or die "Error calling open2: $!";
        print CHLD_IN "$id\n";
        close CHLD_IN;
	my $seqIO = Bio::SeqIO->new(-fh => \*CHLD_OUT, -format => 'fasta');
	my $seq = $seqIO->next_seq();
        close CHLD_OUT;
        wait;
        return $seq;
}

sub cleanup {
	wait;
}

sub get_blast_io {
	my $input_file = shift;
	die "Error lost input file" unless defined $input_file;
	my $blast_db = shift;
	die "Error lost blast database file" unless defined $blast_db;
	my $blast_prog = shift;
	die "Error lost blast program" unless defined $blast_prog;
	my $num_descs = shift;
	die "Error lost number of descriptions" unless defined $num_descs;

	# Call BLAST
	my $pid = open2(\*RDRFH, \*WTRFH, 'blastall',
		'-i',	$input_file,
		'-p',	$blast_prog,
		'-d',	$blast_db,
		'-a',	3,
		'-F',	'F',
		'-v',	$num_descs,
		'-b',	$num_descs
	) or die "Error calling open2: $!";
	close WTRFH or die "Error closing STDIN to blastall process: $!";
#	warn "blastall -i $input_file -p $blast_prog -d $blast_db -a 3 -F F -v $num_descs -b $num_descs";
	
	my $blastio = Bio::SearchIO->new(-fh => \*RDRFH, -format => 'blast') or die "Error creating blastIO: $!";
	return ($blastio,$pid);
}

sub usage {
	die <<EOF 
Usage: $0 -i fasta -d blast_db -p blast_prog -v num_descs
N.B.	blast_prog must be either blastp or blastn
EOF
}

sub ensure_blast_formated {
	my $fasta = shift;
	die "Error lost FASTA for BLAST database" unless defined $fasta;
	my $blast_prog = shift;

	my $formatdb = 'formatdb';

	my @extensions;
	if ($blast_prog eq 'blastp'){
		@extensions = qw( phr pin psq );
	}elsif($blast_prog eq 'blastn'){
		@extensions = qw( nhr nin nsq );
	}else{
		die "Error $blast_prog is an unsupported type of BLAST search";
	}
	my $num_found = 0;
	foreach my $ext (@extensions){
		my $target = join('.',$fasta,$ext);
		$num_found ++ if (-e $target);
	}

	# If there are missing files then format the database
	if ($num_found != @extensions){
		my $protein;
		($blast_prog eq 'blastp') ? $protein = 'T' : $protein = 'F';
		system ($formatdb,'-i',$fasta,'-p',$protein) == 0 or die "Error calling $formatdb -i $fasta -P $protein: $!";
	}
}


