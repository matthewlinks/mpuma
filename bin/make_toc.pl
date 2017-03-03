#!/usr/bin/perl -w
use Getopt::Long;
my ($dir,$output,$VERBOSE,$nCPU,$includeFasta);
my @options = (
	'd=s',	\$dir,
	'o=s',	\$output,
	'v=s',	\$VERBOSE,
	'n=s',	\$nCPU,
	'i=s',	\$includeFasta
);
&GetOptions(@options);

die "Usage: $0 -d dir -o output.toc"
unless (
	defined $dir
	and defined $output

);
$VERBOSE = 0 unless defined $VERBOSE;
$nCPU = 1 unless defined $nCPU;

my %includeIDs;
if (defined $includeFasta){
	open(FILE,"<$includeFasta") or die "Error opening $includeFasta: $!";
	while(<FILE>){
		chomp $_;
		next unless ($_ =~ /^>/);
		my ($id) = split(/\s+/,$_);
		if ($id =~ /^\>(.*)$/){
			$includeIDs{$1}++;
		}
	}
	close FILE or die "Error closing $includeFasta: $!";
}

# setup the output
open(OUTPUT,">$output") or die "Error opening $output for writting: $!";

# parse all the isotig.ace files
        
my $num = 0;
 
my $files = find_files ($dir);
my %toc;
if ((require Parallel::Loops)and($nCPU > 1)){
	my $parallel = Parallel::Loops->new($nCPU);
	$parallel->share(\%toc);
	$parallel->foreach($files,sub {
		my $file = $_;
		my $aceFile = join('/',$dir,$file);
		my $out = $aceFile . '.reads';
		warn "Reading data from $aceFile" if $VERBOSE;
		open(FILE,"<$aceFile") or die "Error opening $aceFile: $!";
		open(OUT,">$out") or die "Error opening $out: $!";
		$toc{$file} = $out;
		while(<FILE>){
			chomp $_;
			if ($_ =~ /^RD/){
				my @parts = split(/\s+/,$_);
				my ($read) = split(/\./,$parts[1]);
				print OUT $read,"\n";
#				push(@{$toc{$file}},$read);
			}else{
				if($_ =~ /^RD/){
					die "Error parsing Read line [$_] from $aceFile";
				}else{
					next;
				}
			}
		}
		close OUT or die "Error closing $out: $!";
		close(FILE) or die "Error closing $aceFile: $!";
	});
}else{
	foreach my $file (@{$files}){
		my $aceFile = join('/',$dir,$file);
		my $out = $aceFile . '.reads';
		$num++;
		warn "$num: Reading data from $aceFile" if $VERBOSE;
		open(FILE,"<$aceFile") or die "Error opening $aceFile: $!";
		open(OUT,">$out") or die "Error opening $out: $!";
		$toc{$file} = $out;
		while(<FILE>){
			chomp $_;
			if ($_ =~ /^RD/){
				my @parts = split(/\s+/,$_);
				my ($read) = split(/\./,$parts[1]);
#				push(@{$toc{$file}},$read);
				print OUT $read,"\n";
			}else{
				if($_ =~ /^RD/){
					die "Error parsing Read line [$_] from $aceFile";
				}else{
					next;
				}
			}
		}
		close OUT or die "Error closing $out: $!";
		close(FILE) or die "Error closing $aceFile: $!";
	}
}

# write out the results
open(OUTPUT,">$output") or die "Error opening $output for writting: $!";
foreach my $f (sort {return $a cmp $b} keys %toc){
	if ($f =~ /^(.*)\.ace$/){
		my $OTU = $1;
		if (defined $includeFasta){
			next unless defined $includeIDs{$OTU};
		}
		my $tocFile = $toc{$f};
		warn "Reading IDs from $tocFile" if $VERBOSE;
		open(FILE,"<$tocFile") or die "Error opening $tocFile: $!";
		while(<FILE>){
			chomp $_;
			print OUTPUT join(' ',$OTU,$_),"\n";
		}
		close FILE or die "Error closing $tocFile: $!";
	}else{
		warn "WARNING: ignoring data for $f because it didn't match isotig00000.ace";
	}
}
close OUTPUT or die "Error closing $output: $!";

exit 0;

sub find_files {
	my $dir = shift;
	die "Error lost dir" unless defined $dir;
	opendir(DIR,$dir) or die "Error opening $dir for reading: $!";
	my @files;
	while(my $file = readdir(DIR)){
#		next unless ($file =~ /^(isotig\d+)\.ace$/);
		next unless ($file =~ /\.ace$/);
		push(@files,$file);
	}
	closedir(DIR) or die "Error closing $dir: $!";
	@files = sort {return $a cmp $b} @files;
	return \@files;
}


