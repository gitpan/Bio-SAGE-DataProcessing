# $Id: test.pl,v 1.1.1.1 2004/05/24 02:52:55 scottz Exp $

use lib './lib';
use Bio::SAGE::DataProcessing;
use strict;

print "Creating new Bio::SAGE::DataProcessing instance...\n";
my $sage = Bio::SAGE::DataProcessing->new();
print "  ...complete.\n";

print "Inputting fake sequence read and quality scores...\n";
my @reads = ( 'AGAGCCATGAAACCGTATGTAGAGAGGGACACATGTAGACAGATAGACACAGATACCATGTATCGCGGCATTACGATCTAGAACATGACGACTATATCGATAGCTAACCATGTTACAAGAGATA' );
my @scores;
for( my $i = 0; $i < length( $reads[0] ); $i++ ) {
  push( @scores, 10 + int(rand( 40 )) );
}
$sage->add_sequence_read( $reads[0] );
$sage->add_quality_scores( \@scores );
$sage->set_quality_cutoff( 0.99 );
print "  ...complete.\n";

print "SAMPLE_DATA:\n";
for( my $i = 0; $i < scalar( @scores ); $i += 30 ) {
  my $j = $i+30;
#  if( $j > scalar( @scores ) ) { $j = scalar( @scores ); }
  if( $j >= scalar( @scores ) ) { $j = scalar( @scores ) - 1; }
  my $subSeq = substr( $reads[0], $i, ( $j-$i+1 ) );
  my @bps = split( //, $subSeq );
  print join( "  ", @bps ) . "\n";
  print join( " ", map { sprintf( "%2i", $_ ) } @scores[$i..$j] ) . "\n";
}

my ( $pDitags, $pDitagScores ) = $sage->extract_ditags();
print "Outputting ditags...\n";
for( my $i = 0; $i < scalar( @$pDitags ); $i++ ) {
  print $$pDitags[$i] . " [" . join( " ", @{$$pDitagScores[$i]} ) . "]\n";
}
print "  ...complete.\n";

my ( $pTags, $pScores ) = $sage->extract_tags();

print "Outputting tags...\n";
for( my $i = 0; $i < scalar( @$pTags ); $i++ ) {
  print $$pTags[$i] . " [" . join( " ", @{$$pScores[$i]} ) . "]\n";
}
#foreach my $tag ( @tags ) { print "    $tag\n"; }
print "  ...complete.\n";

print "Creating a fake fasta sequence file...\n";

# create a fake FASTA sequence file
open( SCRIPT, $0 ) || die( "Couldn't open this script: $0" );
open( FASTA, '>fakedata.fasta' ) || die( "Could create file fakedata.fasta" );
my $bParsing = 0;
while( my $line = <SCRIPT> ) {
  chomp( $line ); $line =~ s/\r//g;
  if( $line =~ /^\/\/ START_READS$/ ) {
    $bParsing = 1; next;
  }
  if( $line =~ /^\/\/ END_READS$/ ) {
    $bParsing = 0; last;
  }
  if( $bParsing == 1 ) {
    print FASTA $line . "\n";
  }
}
close( FASTA );
close( SCRIPT );

print "  ...complete.\n";
print "Inputting fake fasta file...\n";

# read the fake file in
$sage->clear();
$sage->set_quality_cutoff( 0.0 );
open( SEQS, 'fakedata.fasta' ) || die( "Could not open file fakedata.fasta" );
$sage->load_fasta_sequences( *SEQS );
close( SEQS );
unlink( 'fakedata.fasta' );

print "  ...complete (deleted fake sequence file).\n";

print "Outputting top 25 tags...\n";

my %counts = %{$sage->extract_tagcounts()};
my $i = 0;
foreach my $tag ( sort { $counts{$b} <=> $counts{$a} } keys %counts ) {
  $i++;
  if( $i > 10 ) { last; }
  print join( "  ", "  ", $tag, $counts{$tag} ) . "\n";
}
print "    ...etc.\n";

print "  ...complete.\n";

print "Generating ditag length distribution...\n";
$sage->print_ditag_length_distribution();
print "  ...complete.\n";

print "Generating ditag base distribution...\n";
$sage->print_ditag_base_distribution();
print "  ...complete.\n";

print "Generating extra base calculation for TACCTGCAGA...\n";

$sage->print_extra_base_calculation( $sage->get_extra_base_calculation( "TACCTGCAGA" ) );

print "  ...complete.\n";

print "Tests completed successfully.\n";

exit( 0 );

=pod

=for comment

// START_READS
>16971
cattgggcctctagatgatgatacgcttaagacttgtgctttcatgaaagcacaagtttt
gaaaaagtcatgaaagcacaagtgcaataagagcacatgaactaaaaaaaaaattacaga
gcatgaaaaataaaatgccgcaggaacgcatggatacagtaagcctcaacaccacatgta
aaatgtatcagccaggagcccatgttgtatcagaagtgaccaccagtcatgagatgataa
aaacgagagcgcttcatggtggccacggcttaaaagcgtgcatgaagacagtggccagta
aataggcatgttgtaatcgtgattctgcaggtacatggaagcacaagtgttttctgcctc
catgagggcttccaataaacattttacatgtgtgttgtcatattctgcaggtacatgtta
aaaagaaacttttaagctacatgaacccgggaagttttttggagtcatgccctgatttta
ttcattgcagacatggtggccacggccacttgtgctttcatggccaggagctaattttgt
attgcatgagcagatcagcatagacccagtcatgttgttattgccttccttcaccaacat
gttgtaaaaggatcagcaggggccatggtggccacggctagatacaaatcatgctcgagc
ggccgccagtgtgatggatatctgcagaatccagcacactggcggccgttactagtggat
ccgagctcggtacaaag
>16971
cattgggcctctagatgcatgacttgttcgctgacatctttctcatgtttcctgctctgg
tcagcaggggccatggaattttataaggccgtggccaccatgttcatacacctagtttcc
cacagcatggagccttggtgctgggttttgccatgccacaggagaacaaactagaaacat
gtacaagacccggccgtggccaccatggtggccacggcctttgaaaaagtcatgtacctg
cagaaacaaaagcaaacatggcaaaaccctggcagaggaccaacatgtgcatctggtgta
taaacctagccatggaacacgtgccccctcggtttccatgatccttgctgactcgtgtcg
tccatgtggtttattaacacttgtgctttcatgccctgattttagtcagcaggggccatg
agaattttctagggcccgtggtgcatgtaatggtaactattctgcaggtacatgacaact
aaaggcagcccctttccatgcccttgaggagcctccactatcccatgcggntactgtgca
taaggggcctcatgtgtgatcagacttctcctgtggcctgttcattatagggtaaattgg
tttaattggaagcaggaccaacattccatccatgatccttgctgataccaagaaatcatg
ggaagaaaaggttgaacccctgtcatggggaaaaccccaaccattaaaaaatgctccaac
cgcccccccattgtgatgggaaatccgcacaaattcccaa
>16971
cattgggcctctagatgatgcaggagtgtgatgtgatagaatcttgtgtctctatccccg
agaactttccctgtgagcacacgtgagaggacaggcgtgggctttccctgtgagatgtga
gccatggggacagctcctcgtgcacatccatgaaaaaaaaaaaggtgttggcaacatgta
cctgcagaatttgagaggaaacatgtacctgcagaatcacaaaaaaagccatgaagaacg
atcaagagcaggaaacatgtacctgcagaatgctggcccctcgcatgaagccccgagaat
tctgcaggtacatgcagctatttcaagggtgttagtcatgcgcttttgtagttcctctac
aaacatgaagaaagttctgcacccccaacaccatgtgtgggaaatctagacaacctacat
ggtgcctgagagggtacaaagtttacatgtaccgacaacagtgggtatggagacatggag
ggccggtgtttttatttttcatgggcgatcctcgatcagcaaggatcatgatgggatttt
tttctgcacgtacatgtctttacctgaatcttttaagaacatgaggctggggtccacttg
tgctttcatgtccgcgtatataccattatacagtcatgacggtccacgaggtataacttc
acatgttggcctggctgtgaacacaaccatgggttccaaaatgctttatttcatgctcta
gcggc
>16971
ctgtgaatagagtcattttttccgcattaaatcaagttacctacagtcaccacttgagaa
aaactaaagacgggtagccgatagaggccctcctcgcatggcatgggccccccccatggg
atttgagacaccagtgtagcggacatagaaaacgcactctggccatggatggagaaaaga
taaaggagtcgactagtgatcatccctaaaagcttatagtgtcgtcacgccgcgtgctat
acctaccccggcatcgccgcgtcttccccgcaactgccccggctgagcccgatgctgatc
tgtccaccccaggcttctctctaacggaggctggtaagcccgcatcaaccttcctcggag
taacaggtataattatagggaacgtacaatcactaatttccacagaccacgccgggtgca
ccatggatattccgccgtcctctcttacctcgatcaccgggccccccctccaaccacgca
accctacgtttagcgctacagtacaagcacggaagttcacctcgcacgtacccgatacaa
cggcatgctccagccccacatatagcaattaccagacgcgctgtggcccgtatccgccgg
gatccctactacctccggccgcgtttaccccccttcgtagtacagaagatattcgacaat
ccgcgttgagacgctacaaagcttgcccccccacacgccattcccgaagattagagtcta
ttcccgctcgcccgacgcggcgcttcaagtctcgataaaacgcaagtaaaaacatacata
taagctgctctcgctctagacaaccccccgataatccaatttccccgtctggccccgttg
taagtatgcggttataacacctcgcccctctatcattccccgtgctctcgcaaataaacg
gccccaccgtaatagatatttcacctcgtgcgtcaacctacagtcggatcgcaaccggaa
aacataaattacaaatatgcgatacaacattgcaccgtttgccccccgccttacatctta
cccacgccgtct
>16971
cttatattacacccaccactgtaattaccaagcaccctaattanaacnnnnnacaaacac
ccacccccggatgggccctctagatgctgacgcggcaatcgggtatggagaatgccacct
ttcccctttttataagatgatacatactgaaaagtgctttccatgcatttgaaagcgatt
ctctgttcatgatcgctttctactcaagggcagccatgatccttgctgattgcttaacta
acatgcaataaactgaacttgtgcttccatgtaaactgtaaatggggcttgcccatgagg
gcttccaagggtggggccacatgtaggttgtctaacatatacgagtcatgataattcttt
gaagtgaaatcacatgactccttcctcccgagaggtgtcatgcaggccccacccacgcaa
ggaaccatgaaagcacaagtgtggaacccgaacatgcttgcttctgtggcccatcagcac
atgatcaaaggtgtcctccgcaaaccatggagcaacccccgctgacttcccatgtcgaat
cacaaaatcccctgaccatcgcagcataacacaacaccacgaaacactgattttctaaca
agccccagctcactacaggacacacgacccaccacgtatacactcaaataaaaaagaccc
tgccaaaaacaaaaaaccaacacccgcaagttaacttccatgttccctcacgaaatcaaa
accgcccttcctagacaaaagcgcccccccaaaaaaaaaaccagccgcaatccacacccc
ccatcacccctcacacctgaaacttcgcaaaacaataccccttaccccatcaaaacggga
cccccccccg
>16971
cattgggcctctgatgatggtggccacggcctttgaaaaagtcatgtcagatctttgagc
caagtcaccatggntgggccacggctacnctggaaaaancatgggaaaaaaaaatccagc
agactcatggtggccacggcccagcaaggatcatgactaacaccctctctcaacacacat
gcatctttttatgttaaccacaaccatgccattgcactccatcaacctctcatgagaaaa
aaaaaaagggtaaagcccatgctcccccaaaagccgtggccaccatgtactcctggcact
tgtgtacaaacatgtgattctaaatacttaatagggacatgatggttaaaggggggcgct
gctccatgtttaaaaattaggtcctgcttccatgcgcttcctctgaaggcagacaacatg
gctacatctcttttttttttgcatgacttatataaaccagggcattacatgattattttt
ctttaaccacaaccatgtattgtgtgtaaagagcaggaaacatgagaaagatgtcttttt
ggaaaccatggctgtggttaaattctgcaggtacatgctgggatgcagctaagccaaacc
catggtggccacggccttttgactgttcatgtgaatgttattaaagttataatgcatgta
cctgcagaagaacacaaatccatgctctcaatatattctgcaggtacatgttgttttaaa
aatcacgtaagcatgtgggcaaagc
>16971
cgattgggccctctagatgcatgaaaaaaaaaattgattgtgcccatggttctcccttct
gtgctgcagccatgtttcctgctcttgttagaagcccatgtacctgcagaatcttattca
agccatgccgaagggtctggctactgagtcatgatcctgagttaggcccacatcctcatg
atctgttgactttttcaccttacatgcaatggagcttcgcgggggggccatgaaagggaa
actgggcgctgctccatgcagcaaaaaaacctctggacaacatgccactcctccagtttt
attaattcatgcttcctgtacaaagagggaatacatgaccttgctgggctggggataaac
acatgcaggccccaccttctgcaggtacatgtttcctgctctttagacaacctacatgca
gtgaggcgagttagacaacctacatgccactgcactcttttttagctgcatgtgaagaga
attattctgcaggtacatgtttctcttcagcaacttcttagcatggagataaatgaaagc
cccagaccatgtttcctgctctgacggagatcaccatgggcttctaacactcttttataa
catgaatttccccttatcatttatctccatgcctcactgacaggggatacaggtcatgac
tcagcccgnagcaggtaagtcatggtgatgggctctgcagtggtaccatgttgggggttt
cctttgaaaaagtcatggtc
>16971
cattgggcctctagatgcatgatccttgctgaggttttgttctcatggtggccacggccc
tgtagcaaaaatggtgaaggcagtaaggacctatatacatgcagctgtcccaggcaggct
ctccatgcttatttgttagcattggcctccatgtcaggaggctatttttttttatcatgt
tctggacccaacttgtgctttcatggtgctgaatggcaaccatcaacaacatggtttcta
tcaatcagcaggtttacatggatcccaactgtgccaatcagccatgctctgatggaaaat
taaccacaaccatgtaaaccaaattatcacgtaaacaccatggagcagcgccctggtttg
tggtcatgaaccagaggtgtttttttttagccatgtgatcgcggctggccgtggccacca
tgatcttctcattgagtgcagtggcatgatcaagggtgtttttttggaaaccatgctgtt
ggtgattcagcaaggatcatggtggccacggctttttttttgtccatgagcggtggtggg
gagaatgttttcatgaagcagacacttaaggcttcatccatgctgggttaatacccaaga
gatccatgttgtttaaaggttatcaaagcacatgtacctgcagaattgcaatgccagaca
tgtctgtaaaggatggggttcaccatgtggtttgagcacacaagcagagccatgtacctg
cagaatgtccgtgggcatattgctcaaccggcgccaa
>16971
tccttccccccccctcccgcccaccctccccacacccaccnnnnnnnnnnnccaagcggc
ctctaatgctgctcgattcagcttctacatgctcctgcaggcattataaaattccatgta
aaccaaattattaacatttatcatgatccttgctgaataacgcaaggcatgatccttgct
gataggataaatccctgccttttgggctttgaaatagctgcatgctgaaaacacaagaca
tagtttacatgcaagcacaagttccttgatccccatgtgtgaaaaaaagccccagacatc
catggcctcctttcactttttgaaaaccatgcctggccctcattttcctagatcctgcct
ccttgcctgtactggagaatcatgcaaacatcaggcacttgtgcttccatgctaaaagcc
tgtgatccacggtgcatgagaaagatgcctaaacacttcacatgtttcctgctctatatg
acacctgcatgccgaaacactgtgtgccgcattacatgggccttctaacacaatcaatcc
ttattggccaaaccgggcacagatgctgccctccatgctccaacggaccgcccaaatgtc
attcaaacccccccctaccatcctcacactgcacctcctcttatacagcacaaaaatctc
gcccaacccctatatcgttaagttcccctcctcccaacataagacaaacaaaattacctt
cgcgcctcaaccatggtaacctccccccctttccttgacggcaaaaccccccactccgcc
caaataccttcctacccgcaccaccctgcggcgtgccactattttaaccccgtacacctc
t
>16971
cccccccnnnnnnnnnnnnnnnnnnnnnnnnncccnnnnncctnttncnnnnnnnnnnnn
nnnnnnnnnncnnnnnnnccccatggcccttgatcagaccgtttctccgatatnggggtn
nnnnccaaannnnacaattgggcctctaaagcatgctcctgggcgcctccccgnaaactg
aaagcacaagtgcagggagagctcatgactttttcaaaagggtggagtggtgtcatgnga
gataaatgnaacttgtggctttcatggggctctgggcactanccgtttgtgcatggaaac
aagggtacactttattgccatgaccaaggagggatccctcattaacatgttctccataac
ctttccgaaccctcatgacttttttcaaaactgggcttccaaggggggcccacgggattt
ttccttggccatgaaacaactgggaatggggttccccccaggggaaggggtttatccacc
cccggaaaccttggcactttttttccccgccccggggtccccgggggagggggccccaaa
ttttattggccaaaaaaaaaaaaagggggccctttttttccccccccccacacggtgggg
gagggtattttccccggggggggaccccacttttttgggggggcccccctccg
>16971
cattgggcctctagatgcatgtacctgcagaatgtgtgtgaaaggatgtcctgctgccgg
accctcccttcatgctcatatgttagattagtagaagcatgctaagacttcattttatgt
ctatcatggacctgacccagacatctttctcatgctgttggtgatactacacacaagcat
gtaagtagcaaatcagtaaggatcatgtaggttgtctaaaaggccaaatcccatggacaa
aaaaaagactacaggcacatgtaacaaaaatgccctccaggtcatggttaaattctagac
atctttctcatgaaagcacaagtgaacacaggtcagcatgttcacagtggctgaaaattg
ggcatgcttataactatgaaaagaaatgcatgacatcatcgatggagaaagggaacatgg
atgtgcacgaatgtttattttcatgaaaaaaaaaaagggcgctgctccatgtggaatgct
ggtaggtgtatgaacatgtacctgcagaatttcttcatcaacatgattccttcaatcagt
gaaatgccatgagcagatcaggtttgaaaaagtcatgaaagcacaagtgaaagtctaaaa
catgagcttctaccactcangataaagcatggtctgcacctgtggtagaagctcatgctt
ccttngctggnggaaaaccaacatgagaaagatgtcatcctgcaggtacatgttctcctg
ctcctatttcatatacctgctccaaccggccgccaatgttga
>16971
ncattgggcctctagatgcatgtactttctcctttgaatagctgcatgtatctgctgaaa
cccagcaaggtcatggtggccacggctcgtgaggctagcatggagcagcgcccctgcagg
tatgccatgngagaaatatattacaaacacaccatgattatccaggggcaaaacatcacc
atgtagtcaaaagaattctgcaggtacatgaaggaagatggttcaaatatatacatgtag
tggagaaatcagcaggggccatgcacaaacggtaaagcggagccaacatgtacctgcaga
aaacaatggtctgcatgccccatactacagcctgaccaacatgttgtaatcgtgcccaag
tgctccatgcattcaaattctcagcaggggccatgaatcctgtggatttttttttaccat
gtacctgcagaatacttttttccacatgcttccttgcctggtagaaagcgatcatgtacc
tgcagaagctccggcgcaccatgtcatataagtagggcgctgctccatggcctcttcccc
catatacgagtcatgttcaatttgtgtaaagtgagcccatgatgtataataagggcgctg
ctccatgaaagcacaagtgcggctctcggtcatgtcagccttctgtgcagctcctccatg
gatctcattttctgtcagtgaggcatgctcattaagagttcctgtacaccatgcgcctgt
aatcccagaggaccaacatgtacctgca
>16971
ccattgggcctctagatgcatgcccccgtgaagttctgcaggtacatgaatatgctttaa
cttgtgcttccatgactaaaaaggaggtgcgggaggcatgcagctatttcaaggccgtgg
ccaccatgtgttctggagaaactcttcacatcatgtttggacccagccctttaagaaaca
tgcgattcaaaaggtctgttctggcatgaggaaagctgcatatttcaattacatggtaat
agaagtaggtcagctcatcatggccatagtatatcgagacacaccatgctaagacttcaa
agtaacaatgcatgtccagaataaaataagccatcccatggaagcacaagtttctgcagg
tacatggtgtggtggtgtttttgaagggcatgagccactgtgcttaagaagtgagcatga
taggtagaggtgaagaggcgggcatgtacctgcagaaggcccccagcccatggtggccac
ggcccaaatgcacttcatggtggccacggccttgagaggaaacatgtacctttaatagta
catagaggcatgaaagcggggctctctccagaacacatgctgtcaccctgagtttaccta
aacatgtacctgcagaattagacaacctacatggtggccacggctctttgaagaccatgc
acgcaatgctatttttttttttcatggttgtggttaatgaaatagctgcatggaggttag
attcacctttattacatgtattaatcaaa
// END_READS

=cut
