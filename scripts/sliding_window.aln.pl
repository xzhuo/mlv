use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;

my $seqio = Bio::SeqIO->new( '-format' => 'Fasta' , -file => 'MLV.full_MLVs.aln.fa');
my $primary_seq;
my @seq_ary = ();
while((my $seq = $seqio->next_seq())) {
	if ($seq->display_id eq "MEXF2090_chr3") {
		$primary_seq = $seq;
	} else {
		push(@seq_ary, $seq);
	}
}

use Bio::Factory::EMBOSS;
my $factory = Bio::Factory::EMBOSS -> new();
my $needle = $factory->program('needle');
my $needleoutfile = 'out.needle';

$needle->run({-asequence => $primary_seq,
             -bsequence => \@seq_ary,
             -outfile   => $needleoutfile});

my $alns = Bio::AlignIO->new(-format => 'emboss',
                              -file   => $needleoutfile);

while ( my $aln = $alns->next_aln ) {
	my $seq = $aln->get_seq_by_pos(2);
	my $id = $seq->display_id;
	# my $len = $aln->length();
	my $primary_obj = $aln->get_seq_by_id("MEXF2090_chr3");
	my $len = $primary_seq->length();
	for (my $i = 1; $i < $len; $i+=100) {
		# my $aln2 = $aln->slice($i,$i+99);
		my $start = $aln->column_from_residue_number( "MEXF2090_chr3", $i);
		my $last_col = $start+99 < $len ? $i+99 : $len;
		my $end = $aln->column_from_residue_number( "MEXF2090_chr3", $last_col - 10);
		my $aln2 = $aln->slice($start,$end);
		my $perc = $aln2->percentage_identity;
		print "$id\t$i\t$perc\n";
	}
}

# use Bio::Tools::Run::StandAloneBlast;
# my $factory = Bio::Tools::Run::StandAloneBlast->new(-program => 'blastn',
#                                                  -outfile => 'bl2seq.out');

# foreach my $seq(@seq_ary){
# 	my $id = $seq->display_id;
# 	my $bl2seq_report = $factory->bl2seq($primary_seq, $seq);
# 	my $str = Bio::AlignIO->new(-file=> 'bl2seq.out',-format => 'bl2seq');
# 	my $aln = $str->next_aln();
# 	my $len = $aln->length();
# 	for (my $i = 1; $i < $len; $i+=100) {
# 		my $aln2 = $aln->slice($i,$i+99);
# 		my $perc = $aln2->percentage_identity;
# 		print "$id\t$i\t$perc\n"
# 	}
# }