#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $count = 25000;

my $outio = Bio::SeqIO->new(-Fh => \*STDOUT, -format => 'fasta');
foreach my $file (@ARGV){
	my $inio = Bio::SeqIO->new(-file => $file, -format => 'fasta');
	while(my $seq = $inio->next_seq()){
		if ($seq->seq() =~ /^[agct]+$/i){
			# no ambiguity characters
			$outio->write_seq($seq);
		}else{
			my $seqs = get_enumerated_seqs($seq);
			my $total = @{$seqs};
			my $num = 0; 
			foreach my $s (@{$seqs}){
				$num++;
				my $new_seq = Bio::PrimarySeq->new(	-id	=> join('_',$seq->display_id(),$num),
									-seq	=> $s);
				$outio->write_seq($new_seq);
				if (($num % $count) == 0){
					warn "wrote $num / $total";
				}
			}
		}
	}
}

exit 0;

sub get_enumerated_seqs {
	my $seq = shift;
	die "Error lost sequence" unless defined $seq;
	my $id = $seq->display_id();

	my @nucs = split(//,$seq->seq());
	my @sequences;
	push(@sequences,''); #dummy first sequence...
	for(my $i = 0; $i < @nucs; $i++){
		my $numberOfSeqs = @sequences;
		warn "Working on position $i in sequence $id: $numberOfSeqs";
		my $nuc = $nucs[$i];
		die "Error have no sequences" if (@sequences < 1);
		if ($nuc =~ /^[agct]$/i){
			# not ambiguous position so just add this onto the end of each sequence we currently have
			for(my $j = 0; $j < @sequences; $j++){
				$sequences[$j] = $sequences[$j] . $nuc;
			}
		}elsif($nuc =~ /^[xni]$/i){
			# any possible base
			my @pos = qw(a g c t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[m]$/i){
			# any possible base
			my @pos = qw(a c);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[r]$/i){
			# any possible base
			my @pos = qw(a g);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[w]$/i){
			# any possible base
			my @pos = qw(a t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[s]$/i){
			# any possible base
			my @pos = qw(c g);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[y]$/i){
			# any possible base
			my @pos = qw(c t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[k]$/i){
			# any possible base
			my @pos = qw(g t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[v]$/i){
			# any possible base
			my @pos = qw(a c g);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[h]$/i){
			# any possible base
			my @pos = qw(a c t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[d]$/i){
			# any possible base
			my @pos = qw(a g t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}elsif($nuc =~ /^[b]$/i){
			# any possible base
			my @pos = qw(c g t);
			@sequences = @{enum_pos(\@sequences,\@pos)};
		}else{
			die "Error no support for ambiguity == $nuc";
		}
		foreach my $sub (@sequences){
			if ((length $sub) != ($i +1)){
				die "Error problem $sub is not $i +1 bps long";
			}
		}
	}
	return \@sequences;
}

sub enum_pos {
	my $seqs = shift;
	die "Error lost sequences" unless defined $seqs;
	my $pos = shift;
	die "Error lost posibilities" unless defined $pos;
	
	# copy the existing sequences
	my @tmp_sequences = @{$seqs};
	# blow away what existed before
	undef($seqs);

	my @sequences; # empty set to recieve the new sequences
	
	foreach my $n (@{$pos}){
		# create a copy of each sequence for as many nucs as there 
		foreach my $s (@tmp_sequences){
			push(@sequences,$s . $n);
		}
	}
	return \@sequences;
}

				
		

