# *%) $Id: DataProcessing.pm,v 1.1.1.1 2004/05/24 02:52:55 scottz Exp $
#
# Copyright (c) 2004 Scott Zuyderduyn <scottz@bccrc.ca>.
# All rights reserved. This program is free software; you
# can redistribute it and/or modify it under the same
# terms as Perl itself.

package Bio::SAGE::DataProcessing;

=pod

=head1 NAME

Bio::SAGE::DataProcessing - A Perl module 
for processing serial analysis of gene expression 
(SAGE) libraries.

=head1 SYNOPSIS

  use Bio::SAGE::DataProcessing;
  $sage = Bio::SAGE::DataProcessing->new();

  # load sequence reads
  open( READS, "library.fasta" );
  $sage->load_fasta_sequences( *READS );
  close( READS );

  # load Phred quality scores
  open( QUALITY, "library.qual.fasta" );
  $sage->load_fasta_scores( *READS );
  close( QUALITY );

  # set quality filter for 95%
  $sage->set_quality_cutoff( 0.95 );

  # output tags in descending order of expression
  my %tags = %{$sage->extract_tagcounts()};
  open( TAGS, ">library.tags" );
  map { print TAGS join( "\t", $_, $tags{$_} ) . "\n" } sort { $tags{$b} <=> $tags{$a} } keys %tags;
  close( TAGS );

  # tag AAACCGGGTT matches two different genes
  # so 15th base pair may help resolve this
  $sage->print_extra_base_calculation( $sage->get_extract_base_calculation( "AAACCGGGTT" ) );

=head1 DESCRIPTION

This module provides several tools for processing and
analyzing serial analysis of gene expression (SAGE)
libraries.

=head1 README

B<BACKGROUND>

Serial analysis of gene expression (SAGE) is a molecular 
technique for generating a near-global snapshot of a 
cell population’s transcriptome.  Briefly, the technique 
extracts short sequences at defined positions of 
transcribed mRNA.  These short sequences are then paired 
to form ditags.  The ditags are concatamerized to form 
long sequences that are then cloned.  The cloned DNA is 
then sequenced.  Bioinformatic techniques are then 
employed to determine the original short tag sequences, 
and to derive their progenitor mRNA.  The number of times
a particular tag is observed can be used to quantitate
the amount of a particular transcript.  The original 
technique was described by Velculescu et al. (1995) and 
utilized an ~14bp sequence tag.  A modified protocol 
was introduced by Saha et al. (2002) that produced ~21bp 
tags.

B<PURPOSE>

This module facilitates the processing of SAGE data.  
Specifically:

  1. extracting ditags from raw sequence reads.
  2. extracting tags from ditags, with the option to
     exclude tags if the Phred scores (described by
     Ewing and Green, 1998a and Ewing et al., 1998b)
     do not meet a minimum cutoff value.
  3. calculating descriptive values
  4. statistical analysis to determine, where possible, 
     additional nucleotides to extend the length of the 
     SAGE tag (thus facilitating more accurate tag to 
     gene mapping).

Both regular SAGE (14mer tag) and LongSAGE (21mer tag) 
are supported by this module.  Future protocols should
be configurable with this module.

B<REFERENCES>

  Velculescu V, Zhang L, Vogelstein B, Kinzler KW. (1995) 
  Serial analysis of gene expression. Science. 270:484-487.

  Saha S, Sparks AB, Rago C, Akmaev V, Wang CJ, Vogelstein B, 
  Kinzler KW, Velculescu V. (2002) Using the transcriptome 
  to annotate the genome. Nat. Biotechnol. 20:508-512.

  Ewing B, Hillier L, Wendl MC, Green P. (1998a) Base-calling
  of automated sequencer traces using phred. I. Accuracy
  assessment. Genome Res. 8:175-185.

  Ewing B, Green P. (1998b) Base-calling of automated
  sequencer traces using phred. II. Error probabilities.
  Genome Res. 8:186-194.

=head1 INSTALLATION

Follow the usual steps for installing any Perl module:

  perl Makefile.PL
  make test
  make install

=head1 PREREQUISITES

This module requires the C<Statistics::Distributions> package.

=head1 CHANGES

  1.00 2004.05.23 - Initial release.

=cut

use strict;
use diagnostics;
use vars qw( $VERSION @ISA @EXPORT @EXPORT_OK $PROTOCOL_SAGE $PROTOCOL_LONGSAGE $DEBUG $ENZYME_NLAIII $ENZYME_SAU3A );

require Exporter;
require AutoLoader;

@ISA = qw( Exporter AutoLoader );
@EXPORT = qw();
$VERSION = "1.00";

use Statistics::Distributions;

my $PACKAGE = "Bio::SAGE::DataProcessing";

=pod

=head1 VARIABLES

B<Globals>

=over 2

I<$PROTOCOL_SAGE>

  Hashref containing protocol parameters for the
  regular/original SAGE protocol (see set_protocol
  documentation for more information).

I<$PROTOCOL_LONGSAGE>

  Hashref containing protocol parameters for the
  LongSAGE protocol (see set_protocol documentation 
  for more information).

I<$ENZYME_NLAIII> = 'CATG'

  Constant denoting the recognition sequence for NlaIII.

I<$ENZYME_SAU3A> = 'GATC'

  Constant denoting the recognition sequence for Sau3a.

=back

B<Settings>

=over 2

I<$DEBUG = 0>

  Prints debugging output if value if >= 1.

=back

=cut

my @ignoreTags = ( 'TCCCTATTAA', 'TCCCCGTACA' ); # linker derived sequences
my @ignoreLongTags = ( 'TCGGACGTACATCGTTA', 'TCGGATATTAAGCCTAG' ); # linker derived sequences

my %params_sage = ( 'MINIMUM_DITAG_LENGTH' => 29,
                    'MAXIMUM_DITAG_LENGTH' => 32,
                    'TAG_LENGTH'           => 10,
                    'IGNORE_TAGS'          => \@ignoreTags );
my %params_longsage = ( 'MINIMUM_DITAG_LENGTH' => 40,
                        'MAXIMUM_DITAG_LENGTH' => 46,
                        'TAG_LENGTH'           => 17,
                        'IGNORE_TAGS'          => \@ignoreLongTags );

$PROTOCOL_SAGE = \%params_sage;         # regular SAGE (14-mer tags)
$PROTOCOL_LONGSAGE = \%params_longsage; # LongSAGE (21-mer tags)

$ENZYME_NLAIII = "CATG";
$ENZYME_SAU3A  = "GATC";

$DEBUG = 0; # set this flag to non-zero to enable debugging messages

=pod

=head1 CLASS METHODS

=cut

#######################################################
sub new {
#######################################################
=pod

=head2 new [$protocol], [$enzyme]

Constructor for a new Bio::SAGE::DataProcessing object.
By default, caching is enabled.

B<Arguments>

I<$protocol> (optional)

  Either $PROTOCOL_SAGE or $PROTOCOL_LONGSAGE. Other
  values will result in an error. The default is
  $PROTOCOL_SAGE.

I<$enzyme> (optional)

  The anchoring enzyme recognition site. This is
  typically NlaIII (CATG) or Sau3a (GATC). The 
  default is the recognition sequence "CATG" (NlaIII).

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new( $Bio::SAGE::DataProcessing::PROTOCOL_LONGSAGE, "CATG" );

=cut

    my $this = shift;
    my $protocol = shift || $PROTOCOL_SAGE;
    my $enzyme = shift || "CATG";
    my $class = ref( $this ) || $this;
    my $self = {};
    bless( $self, $class );

    # set instance variables
    $self->{'enzyme'} = $enzyme;
    $self->set_quality_cutoff( 0 );
    $self->set_protocol( $protocol );
    $self->{'enable_cache'} = 1;

    # prepare cache
    my %cacheHash;
    $self->{'cache'} = \%cacheHash;

    return $self;
    
}

#######################################################
sub load_ditag_sequences {
#######################################################
=pod

=head2 load_ditag_sequences $handle

Creates an array of ditags (compatible with methods
like extract_tags or get_extra_base_calculation) from
an input source.

The expected format of the input is a single ditag
sequence on each line.

NOTE: This method can be called on an instantiated
object, as well as statically.

B<Arguments>

I<$handle>

  A handle to an open file, STDIN, etc.

B<Usage>

  my @ditags = @{Bio::SAGE::DataProcessing::load_ditag_sequences( *STDIN )};

  # this also works
  my $sage = Bio::SAGE::DataProcessing->new();
  my @ditags = @{$sage->load_ditag_sequences( *STDIN )};

=cut

    if( ref( $_[0] ) &&
        ref( $_[0] ) eq $PACKAGE ) {
        shift; # called from instance, ignore reference to object
    }

    my $handle = shift || die( $PACKAGE . "::load_ditag_sequences no data handle provided." );

    my @ditags;
    while( my $line = <$handle> ) {
        chomp( $line ); $line =~ s/\r//g;
        push( @ditags, $line );
    }

    return \@ditags;

}

#######################################################
sub load_ditag_scores {
#######################################################
=pod

=head2 load_ditag_scores $handle

Creates an arrayref to ditag Phred scores 
(compatible with extract_tags) from an input source.

The expected format of the input is a single line of
Phred scores separated by spaces.  Since this array
is usually paired with ditag sequences, the order
of the score file should match the order of the
ditag sequence input source (see load_ditag_sequences).

NOTE: This method can be called on an instantiated
object, as well as statically.

B<Arguments>

I<$handle>

  A handle to an open file, STDIN, etc.

B<Usage>

  my @scores = @{Bio::SAGE::DataProcessing::load_ditag_scores( *STDIN )};

  # this also works
  my $sage = Bio::SAGE::DataProcessing->new();
  my @scores = @{$sage->load_ditag_scores( *STDIN )};

  map { print join( " ", @{$_} ) . "\n" } @scores;

=cut

    if( ref( $_[0] ) &&
        ref( $_[0] ) eq $PACKAGE ) {
        shift; # called from instance, ignore reference to object
    }

    my $handle = shift || die( $PACKAGE . "::load_ditag_scores no data handle provided." );

    my @scores;
    while( my $line = <$handle> ) {
        chomp( $line ); $line =~ s/\r//g;
        my @arr = split( /\s+/, $line );
        push( @scores, \@arr );
    }

    return \@scores;

}

=pod

=head1 INSTANCE METHODS

=cut

#######################################################
sub clear {
#######################################################
=pod

=head2 clear

Clears this object of all added sequences
and quality scores.

B<Arguments>

  None.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  $sage->clear();

=cut

    my $this = shift;

    $this->{'sequences'} = ();
    $this->{'scores'} = ();
    $this->__clearCache();

    return;

}

#######################################################
sub is_caching_enabled {
#######################################################
=pod

=head2 is_caching_enabled

Gets whether this object caches processed ditags, tags,
and other calculations.

B<Arguments>

  None.

B<Returns>

  Returns 1 if caching is enabled, 0 otherwise.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  if( $sage->is_caching_enabled() ) {
    print STDERR "Caching enabled.\n";
  }

=cut

  my $this = shift;

  return $this->{'enable_cache'};

}

#######################################################
sub set_caching_enabled {
#######################################################
=pod

=head2 set_caching_enabled $boolean

Sets whether this object should cache processed ditags,
tags, and the results of other calculations so 
subsequent processing is accelerated.  By default,
new objects have caching enabled.

B<Arguments>

I<$boolean>

  Any value >= 1 is considered TRUE, 0 is considered
  FALSE.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  $sage->set_caching_enabled( 1 );

=cut

    my $this = shift;
    my $boolean = shift;

    if( !defined( $boolean ) ) { die( $PACKAGE . "::set_caching_enabled no argument specified.\n" ); }

    if( $boolean < 0 ) { die( $PACKAGE . "::set_caching_enabled invalid argument: $boolean\n" ); }

    if( $boolean > 0 ) { $boolean = 1; }

    if( $boolean == $this->{'enable_cache'} ) {
        return; # no change
    }
  
    if( $boolean == 0 ) {
       # disable caching
       # - remove any currently cached data
       %{$this->{'cache'}} = ();
    }

    $this->{'enable_cache'} = $boolean;

    return;

}

#######################################################
sub get_enzyme {
#######################################################
=pod

=head2 get_enzyme

Gets the current anchoring enzyme recognition site.

B<Arguments>

  None.

B<Returns>

  The current anchoring enzyme recognition site. By
  default, this will be 'CATG', the NlaIII recognition
  site.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  print "ENZYME_SITE: " . $sage->get_enzyme() . "\n";

=cut
    
    my $this = shift;

    return $this->{'enzyme'};
    
}

#######################################################
sub set_enzyme {
#######################################################
=pod

=head2 set_enzyme $enzyme

Sets the anchoring enzyme recognition site.

B<Arguments>

I<$enzyme>

  The anchoring enzyme recognition site. This is
  typically NlaIII (CATG) or Sau3a (GATC). The 
  default is the recognition sequence "CATG" (NlaIII).

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  $sage->set_enzyme( $Bio::SAGE::DataProcessing::ENZYME_SAU3A );

=cut
    
    my $this = shift;
    my $enzyme = shift || die( $PACKAGE . "::set_enzyme no enzyme recognition site specified." );

    if( $this->{'enzyme'} ne $enzyme ) {
      # enzyme is changing, so clear cache
      $this->__clearCache();
    }

    $this->{'enzyme'} = $enzyme;
    
}

#######################################################
sub get_protocol {
#######################################################
=pod

=head2 get_protocol

Gets the SAGE protocol parameters this object is currently 
set to use.

B<Arguments>

  None.

B<Returns>

  A hashref containing protocol parameters.
  See the set_protocol documentation for more
  information.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my %params = %{$sage->get_protocol()};
  print "CURRENT PARAMETERS:\n";
  map { print join( "\t", $_, $params{$_} ) . "\n" } keys %params;  

=cut
    
    my $this = shift;
    
    return $this->{'protocol'};
    
}

#######################################################
sub set_protocol {
#######################################################
=pod

=head2 set_protocol \%protocolParams

Sets the SAGE protocol parameters that should be used 
by this object to process sequence data.

B<Arguments>

I<\%protocolParams>

  A hashref containing specifics of the protocol. Two 
  pre-made parameter sets are available: $PROTOCOL_SAGE 
  (regular SAGE) and $PROTOCOL_LONGSAGE (LongSAGE).

  The required hash keys:

    MINIMUM_DITAG_LENGTH | The shortest length a valid
                           ditag can be (the length
                           should include the anchoring
                           enzyme site sequences).
    MAXIMUM_DITAG_LENGTH | The longest length a valid
                           ditag can be (the length
                           should include the anchoring
                           enzyme site sequences).
    TAG_LENGTH           | The expected tag length (the
                           length should NOT include the
                           anchoring enzyme site sequence).
    IGNORE_TAGS          | An arrayref listing tag
                           sequences that should be
                           ignored during tag extraction
                           (i.e. linker-derived sequences).

  The parameters for the default configurations are:

                           +----------------+--------------------+
                           | $PROTOCOL_SAGE | $PROTOCOL_LONGSAGE |
    +----------------------+----------------+--------------------+     
    | MINIMUM_DITAG_LENGTH |       29       |         40         |
    +----------------------+----------------+--------------------+     
    | MAXIMUM_DITAG_LENGTH |       32       |         46         |
    +----------------------+----------------+--------------------+     
    | TAG_LENGTH           |       10       |         14         |
    +----------------------+----------------+--------------------+     
    | IGNORE_TAGS          |   TCCCTATTAA   | TCGGACGTACATCGTTA  |
    |                      |   TCCCCGTACA   | TCGGATATTAAGCCTAG  |
    +----------------------+----------------+--------------------+     

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  $sage->set_protocol( $Bio::SAGE::DataProcessing::PROTOCOL_LONGSAGE );

  # alternatively:
  my %params = ( 'MINIMUM_DITAG_LENGTH' => 31,
                 'MAXIMUM_DITAG_LENGTH' => 34,
                 'TAG_LENGTH'           => 11,
                 'IGNORE_TAGS'          => \( 'TCCCTATTAA', 'TCCCCGTACA' ) );
  $sage->set_protocol( \%params );

=cut
    
    my $this = shift;
    my $protocol = shift;
    
    my %params = %{$protocol};

    if( !defined( $params{'MINIMUM_DITAG_LENGTH'} ) ) {
      die( $PACKAGE . "::set_protocol MINIMUM_DITAG_LENGTH not specified." );
    }
    if( !defined( $params{'MAXIMUM_DITAG_LENGTH'} ) ) {
      die( $PACKAGE . "::set_protocol MAXIMUM_DITAG_LENGTH not specified." );
    }
    if( !defined( $params{'TAG_LENGTH'} ) ) {
      die( $PACKAGE . "::set_protocol TAG_LENGTH not specified." );
    }
    
    $this->{'protocol'} = \%params;
    
    $this->__clearCache(); # clear any previous analysis

    return;
    
}

#######################################################
sub get_quality_cutoff {
#######################################################
=pod

=head2 get_quality_cutoff

Gets the quality cutoff for a tag to be considered
valid during extraction.

B<Arguments>

  None.

B<Returns>

  The minimum quality (0.0 <= qual < 1.0) required
  for a tag to be considered valid.  By default,
  this value is 0.0.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  print "QUALITY_CUTOFF: " . $sage->get_quality_cutoff();

=cut
    
    my $this = shift;

    return $this->{'quality_cutoff'};
    
}

#######################################################
sub set_quality_cutoff {
#######################################################
=pod

=head2 set_quality_cutoff $cutoff

Sets the quality cutoff for a tag to be considered
valid during extraction.

B<Arguments>

I<$cutoff>

  A value from 0.00 and <1.00.  For example, a quality 
  cutoff where the cumulative probability of correct 
  sequence should be 99% would set this function with
  the value 0.99.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  $sage->set_quality_cutoff( 0.99 );

=cut
    
    my $this = shift;
    my $cutoff = shift;

    if( !defined( $cutoff ) ) { die( $PACKAGE . "::set_quality_cutoff no cutoff value was passed." ); }

    if( $cutoff < 0.0 || $cutoff >= 1.0 ) {
      die( $PACKAGE . "::seq_quality_cutoff cutoff value $cutoff is not 0 <= \$cutoff < 1.0\n" );
    }
    
    $this->{'quality_cutoff'} = $cutoff;

    $this->__clearCache(); # clear any previous analysis
    
    return;
    
}

#######################################################
sub add_sequence_read {
#######################################################
=pod

=head2 add_sequence_read $sequence

Adds a sequence read to current collection of reads.

B<Arguments>

I<$sequence>

  A valid nucleotide sequence. The only valid characters
  are A,C,G,T,N,X.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  $sage->add_sequence_read( $read );

=cut

    my $this = shift;
    my $sequence = shift;

    if( $sequence !~ /^[ACGTacgtNnXx]+$/ ) {
        die( $PACKAGE . "::add_sequence_read invalid sequence " . $sequence . "\n" );
    }

    # adding new data, so cache will no longer
    # be valid - clear it
    $this->__clearCache();

    push( @{$this->{'sequences'}}, $sequence );  

}

#######################################################
sub add_quality_scores {
#######################################################
=pod

=head2 add_quality_scores \@scores

Adds Phred scores for the last sequence added using
add_sequence_read.

B<Arguments>

I<\@scores>

  A arrayref of Phred scores. This method must be called 
  after each addition of a nucleotide  sequence. The size 
  of the array must match the length of the last added 
  sequence. A valid Phred score ranges from 0-99.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my $sequence = 'ACGTACGTACGT';
  my @scores = ( 20, 30, 34, 32, 23, 25, 14, 45, 12, 34, 23, 25 );
  $sage->add_sequence_read( $read );
  $sage->add_quality_scores( \@scores );

=cut

    my $this = shift;
    my $pQuality = shift;

# TODO: Debug
#    if( scalar( @{$this->{'sequences'}} ) == 0 ) {
#        die( $PACKAGE . "::add_quality_scores no sequence.\n" );
#    }

# TODO: Debug
#    if( scalar( @{$this->{'sequences'}} ) != scalar( @{$this->{'quality'}} ) {
#        die( $PACKAGE . "::add_quality_scores must be called after corresponding sequence is added.\n" );
#    }

# TODO: Debug
#    my $seqLength = length( ${$this->{'sequences'}}[scalar(@{$this->{'sequences'}})-1] );
#    if( $seqLength != scalar( @$pQuality ) ) {
#        die( $PACKAGE . "::add_quality_scores quality scores and sequence length don't match.\n" );
#    }

    map { 
        if( $_ < 0 || $_ > 99 ) { die( $PACKAGE . "::add_quality_scores one or more Phred scores was not between 0-99.\n" ); } 
    } @$pQuality;

    # adding new data, so cache will no longer
    # be valid - clear it
    $this->__clearCache();

    push( @{$this->{'scores'}}, $pQuality );

}

#######################################################
sub load_fasta_sequences {
#######################################################
=pod

=head2 load_fasta_sequences $handle

Loads sequence data in FASTA format.

An example of FASTA format sequence:

  >SEQUENCE0001
  ACAGATAGACAGAGATATAGAGACATATTTAGAGACAAATCGCGCAGGCGCGCGACATA
  TGACTAGTTTATATCATCAGTATTAGCGATTATGACTATTATATATTACTGATTGATCT
  ATAGCGCGATTATATCTATCTATGCATTCGATCATGCTATTATCGTATACTACTGCTAG
  AGAGGACGACGCAGGCAGCGATTATATCTATTTATA
  >SEQUENCE0002
  CGCGACGCATGTCAGTAGCTAGCTGCGCCCGAATATATATCGTCATACGGATTCGTAGC
  CCCCCGCGGAGTCTGATTATATCTGATT

B<Arguments>

I<$handle>

  A Perl handle to the input data.

B<Returns>

  The number of sequences read.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my $seqFile;
  open( $seqFile, "reads.fasta" );
  print "LOADED " . $sage->load_fasta_sequences( $seqFile ) . " SEQUENCES.\n";
  close( $seqFile );

=cut

  my $this = shift;
  my $fh = shift || die( $PACKAGE . "::load_fasta_sequences No handle to FASTA content given.\n" );

  my $currid = '';
  my $currSeq = '';

  my $nLoaded = 0;

  while( my $line = <$fh> ) {
    chomp( $line ); $line =~ s/\r//g;
    if( $line =~ /^>(.*?)$/ ) {
        my $thisid = $1;
        if( $currid ne '' ) {
            push( @{$this->{'sequences'}}, $currSeq );
            $nLoaded++;
            $currSeq = '';
        }
        $currid = $thisid;
        next;
    }
    if( $line =~ /^[^\d]+$/ ) {
        $currSeq .= $line; next;
    }
    die( $PACKAGE . "::load_fasta_sequences unrecognized line: " . $line . "\n" );
  }
  if( $currid ne '' ) {
    push( @{$this->{'sequences'}}, $currSeq );
    $nLoaded++;
    $currSeq = '';
  }

  return $nLoaded;

}

#######################################################
sub load_fasta_scores {
#######################################################
=pod

=head2 load_fasta_scores $handle

Loads Phred scores in FASTA format.

An example of FASTA format quality data:

  >SEQUENCE0001
  11 17 18 16 19 17 19 19 16 19 19 16 11 10 9 15 10 12 24 24 35 29 29 39 40 40 40 40 37 37 46 46 40 40 40 40 56 56 56 56 35 35 35 35 35 35 56 40 40 46 40 40 39 39 35 39 56 56 51 51
  51 51 51 56 35 35 35 35 35 35 40 40 51 45 51 51 39 39 39 39 39 39 40 40 40 40 40 40 56 56 56 56 56 46 46 40 39 39 39 45 45 45 56 56 56 56 56 56 56 56 40 39 39 39 39 35 39 39 39 39
  45 56 45 45 45 45 51 35 39 39 39 39 39 40 40 40 40 40 51 56 56 40 40 40 40 40 43 56 56 56 43 35 35 35 35 35 43 45 45 45 45 45 45 51 51 51 51 51 51 56 56 56 56 56 56 51 51 51 56 56
  7 7 9 9 11 10 13 11 10 8 10 10 8 8 8 10 10
  >SEQUENCE0002
  12 15 17 17 19 15 15 15 13 19 17 17 12 16 11 19 13 24 24 35 35 35 37 37 39 56 56 56 56 56 51 39 32 29 29 29 29 32 56 56 56 35 35 35 35 35 35 56 56 56 56 56 56 56 56 56 56 56 56 40
  40 40 46 46 40 51 40 40 40 40 40 40 51 39 39 35 35 35 35 40 40 51 45 45 45 45 51 51 56 56 56 56 56 45 45 45 45 51 51 45 45 45 40 40 40 40 40 40 40 40 40 40 56 56 56 56 56 56 51 51
  15 13 19 17 17 12 16 11 19 13 

B<Arguments>

I<$handle>

  A Perl handle to the input data.

B<Returns>

  The number of score sets read.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my $qualFile;
  open( $qualFile, "reads.qual.fasta" );
  print "LOADED " . $sage->load_fasta_scores( $qualFile ) . " QUALITY SETS.\n";
  close( $qualFile );

=cut

    my $this = shift;
    my $fh = shift || die( $PACKAGE . "::load_fasta_scores No handle to FASTA content given.\n" );

    my $currid = '';
    my $currSeq = '';

    my $nLoaded = 0;

    while( my $line = <$fh> ) {
        chomp( $line ); $line =~ s/\r//g;
        if( $line =~ /^>(.*?)$/ ) {
            my $thisid = $1;
            if( $currid ne '' ) {
                my @scores = split( /\s+/, $currSeq );
                push( @{$this->{'scores'}}, \@scores );
                $nLoaded++;
                $currSeq = '';
            }
            $currid = $thisid;
            next;
        }
        if( $line =~ /^[\d\s]+$/ ) {
            if( $currSeq ne '' ) { $currSeq .= ' '; }
            $currSeq .= $line;
            next;
        }
        die( $PACKAGE . "::load_fasta_scores unrecognized line: " . $line . "\n" );
    }
    if( $currid ne '' ) {
        my @scores = split( /\s+/, $currSeq );
        push( @{$this->{'scores'}}, \@scores );
        $nLoaded++;
        $currSeq = '';
    }

    if( $DEBUG >= 1 ) {
        print STDERR "Loaded $nLoaded scores.\n";
    }

    return $nLoaded;

}

#######################################################
sub get_number_sequences {
#######################################################
=pod

=head2 get_number_sequences

Gets the number of sequences that have been added
to this object.

B<Arguments>

  None.

B<Returns>

  The number of sequences that have been added to
  this object with add_sequence_data( $sequence ).

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  foreach my $sequence ( @sequences ) {
    $sage->add_sequence_data( $sequence );
  }
  print "Added " . $sage->get_number_sequences() . " sequences.\n";

=cut

    my $this = shift;

    return scalar( @{$this->{'sequences'}} );

}

#######################################################
sub extract_ditags {
#######################################################
=pod

=head2 extract_ditags [\@sequences], [\@scores]

Extracts valid ditags from sequence reads and returns
an arrayref of the list of ditags. A valid ditag must
be within the minimum and maximum allowable length.
The allowable length is dependent on the currently
selected SAGE protocol.

B<Arguments>

I<\@sequences> (optional)

  An arrayref of sequence reads. If this argument is not 
  supplied, the reads that have been added to this object 
  with add_sequence_read are used NOTE: if caching is 
  turned on, it will only be helpful for ditag extraction 
  operations on reads stored in this object (i.e. when no 
  \@sequences argument is provided).

I<\@scores> (optional)

  An arrayref containing  score arrays corresponding
  to the specified ditag sequences.  This must be
  supplied if get_quality_cutoff > 0.

B<Returns>

  An arrayref containing ditag sequences and, if
  quality data is available and quality filtering is 
  enabled (get_quality_cutoff > 0), a arrayref
  containing arrays of Phred scores. 

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @reads = ( 'ACGTAGACATAGACAAGAGATATAGA',
                'GATAGACAAAGGAAGATTACAAGATTAT' );

  # extract and print ditags
  map { print $_ . "\n" } @{$sage->extract_ditags( \@reads )};

  # alternatively:
  map { $sage->add_sequence_read( $_ ) } @reads;
  map { print $_ . "\n" } @{$sage->extract_ditags()};

  # for quality
  $sage->set_quality_cutoff( 0.99 );
  $sage->add_quality_scores( \( 20, 30, 34, ..., 30 ),
                             \( 34, 23, 12, ..., 21 ) );
  my( $pDitags, $pScores ) = $sage->extract_ditags();
  map { print $_ . "\n" } @{$pDitags}; # print ditags
  map { print join( " ", @{$_} ) . "\n" } @{$pScores}; # print scores

=cut

    my $this = shift;
    my $pSequences = shift || $this->{'sequences'};

    my $bQuality = 0;
    my $pScores;

    # are we using internal data?
    if( $pSequences == $this->{'sequences'} ) {

        # are we prepared to use quality?
        if( $this->get_quality_cutoff() > 0.0 ) {
            $bQuality = 1;
            if( scalar( @{$this->{'sequences'}} ) != scalar( @{$this->{'scores'}} ) ) {
                die( $PACKAGE . "::extract_ditags quality filtering is set but number of Phred score sets does not match number of sequences (did you add Phred scores?)." );
            }
            $pScores = $this->{'scores'};
        }

        # do we have a cache of this data?
        if( $this->is_caching_enabled() && 
            defined( $this->{'cache'}->{'ditags'} ) && 
            scalar( @{$this->{'cache'}->{'ditags'}} ) > 0 ) {
            if( $bQuality == 1 ) {
                return ( \@{$this->{'cache'}->{'ditags'}}, \@{$this->{'cache'}->{'ditag_scores'}} );
            }
            return \@{$this->{'cache'}->{'ditags'}};
        }

    } else {

        if( $this->get_quality_cutoff() > 0.0 ) {
            $bQuality = 1;
            # did user also pass quality scores?
            $pScores = shift || die( $PACKAGE . "::extract_ditags wants to filter with cutoff " .
                                     $this->get_quality_cutoff() . " but no scores were provided." );
            if( scalar( @$pSequences ) != scalar( @$pScores ) ) {
                die( $PACKAGE . "::extract_ditags quality filtering is set but number of Phred score sets does not match number of sequences (did you provide scores as an argument?)." );
            }
        }

    }
    
    my @fragments;
    my $site = $this->{'enzyme'};
    
    my $minLength = $this->{'protocol'}->{'MINIMUM_DITAG_LENGTH'};
    my $maxLength = $this->{'protocol'}->{'MAXIMUM_DITAG_LENGTH'};

    if( $DEBUG >= 1 ) {
        print STDERR "MINLENGTH=$minLength\n";
        print STDERR "MAXLENGTH=$maxLength\n";
    }

    my @scores;

    for( my $x = 0; $x < scalar( @$pSequences ); $x++ ) {
        my $read = $$pSequences[$x];
        $read = uc( $read ); # make bases uppercase
        my $pos = -1;
        my @positions;
        while( ( $pos = index( $read, $site, $pos ) ) > -1 ) {
            push( @positions, $pos );
            $pos++;
        }
        for( my $i = 0; $i < scalar( @positions )-1; $i++ ) {
            my $fragLength = $positions[($i+1)] + length($site) - $positions[$i] + 1;
            if( $minLength > $fragLength || $fragLength > $maxLength ) { next; } # too long or short
            push( @fragments, substr( $read, $positions[$i], $positions[($i+1)]+length($site)-$positions[$i] ) );
            if( $bQuality ) {
                my @arr = @{$$pScores[$x]};
                my @sec = @arr[$positions[$i]..($positions[($i+1)]-1+length($site))];
                push( @scores, \@sec );
            }
        }
    }
    
    if( $pSequences == $this->{'sequences'} && $this->is_caching_enabled() ) {
        # Cache results
        @{$this->{'cache'}->{'ditags'}} = ();
        push( @{$this->{'cache'}->{'ditags'}}, @fragments );
        push( @{$this->{'cache'}->{'ditag_scores'}}, @scores );
    }

    if( $bQuality ) {
        return ( \@fragments, \@scores );
    }

    return \@fragments;
    
}

#######################################################
sub extract_tags {
#######################################################
=pod

=head2 extract_tags [\@ditags], [\@scores]

Extracts valid tags from ditags and returns a reference
to the list of tags. A valid tag contains only
A,C,G,T (extracted tags containing N or X are ignored).
The length of extracted tags is dependent on the
currently selected SAGE protocol.

The right-hand SAGE tag is reverse complemented. For
example:

  CATGAAACCGTATGTAGAGAGGGACACATG

becomes:

  AAACCGTATG
  TGTCCCTCTC

NOTE: Tags that have been set to be ignored (see the
set_protocol method) will not be included in the
output.

B<Arguments>

I<\@ditags> (optional)

  An arrayref containing ditag sequences. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used. (NOTE: if caching is 
  turned on, it will only be applied to tag extraction 
  operations on reads/ditags stored in this object 
  (ie. when no \@ditags argument is provided)).

I<\@scores> (optional)

  An arrayref containing score arrays corresponding
  to the specified ditag sequences.

B<Returns>

  An arrayref containing tag sequences. If quality 
  data is available and quality filtering is enabled 
  (get_quality_cutoff > 0) an arrayref containing 
  Phred score arrays for each tag is also returned.
  (NOTE: The Phred score array will include scores
  for the anchoring enzyme site.  For example, a
  regular SAGE library constructed with the NlaIII
  enzyme would return a 10-mer tag, and a 14 element
  score array that includes the scores for the
  leading CATG nucleotides).

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @reads = ( 'ACGTAGACATAGACAAGAGATATAGA',
                'GATAGACAAAGGAAGATTACAAGATTAT' );

  my @ditags = @{$sage->extract_ditags( \@reads )};

  # extract and print tags
  map { print $_ . "\n" } @{$sage->extract_tags( \@ditags )};

  # alternatively:
  map { $sage->add_sequence_read( $_ ) } @reads;
  map { print $_ . "\n" } @{$sage->extract_tags()};

  # alternatively:
  my @scores = ( [ 20, 30, 34, ..., 30 ],
                 [ 34, 23, 12, ..., 21 ] );
  $sage->set_quality_cutoff( 0.99 );
  my ( $pDitags, $pScores ) = $sage->extract_tags( \@ditags, \@scores );
  for( my $i = 0; $i < scalar( @$pDitags ); $i++ ) {
    print $$pDitags[$i] . "\n";
    print join( " ", $$pScores[$i] ) . "\n";
  }

=cut

    my $this = shift;
    my $pDitags = shift;
    my $pScores = shift;

    if( defined( $pDitags ) && !defined( $pScores ) && $this->get_quality_cutoff() > 0.0 ) {
        die( $PACKAGE . "::extract_tags wants to use quality cutoff but only ditag sequences were provided (need scores too)." );
    }

    my $bInternal = 0;
    my $bQuality = 0;
    # are no ditags provided?
    if( !defined( $pDitags ) ) {

        # are we prepared to use quality?
        if( $this->get_quality_cutoff() > 0.0 ) {
            $bQuality = 1;
            ( $pDitags, $pScores ) = $this->extract_ditags();
            if( scalar( @$pDitags ) != scalar( @$pScores ) ) { # paranoia check
                die( $PACKAGE . "::extract_tags internal error - sequence and score arrays aren't same size." );
            }
        } else {
            $pDitags = $this->extract_ditags();
        }

        $bInternal = 1;
        # do we have a cache of this data?
        if( $this->is_caching_enabled() && 
            defined( $this->{'cache'}->{'tags'} ) &&
            scalar( @{$this->{'cache'}->{'tags'}} ) > 0 ) {
            return \@{$this->{'cache'}->{'tags'}};
        }

    }

    my $tagLength = $this->{'protocol'}->{'TAG_LENGTH'};
    my $enzyme = $this->{'enzyme'};
 
    my @tags;
    my @scores;
    
    for( my $x = 0; $x < scalar( @$pDitags ); $x++ ) {
        my $ditag = $$pDitags[$x];

        my $tag1 = substr( $ditag, length( $enzyme ), $tagLength );
        my $ignoreTag1 = 0;
        foreach my $ignoreTag ( @{$this->{'protocol'}->{'IGNORE_TAGS'}} ) {
            if( $ignoreTag eq $tag1 ) { $ignoreTag1 = 1; next; }
        }
        if( $ignoreTag1 == 0 ) {
            if( $tag1 =~ /^[ACGT]+$/ ) { # must contain only A,C,G,T
                push( @tags, substr( $ditag, length( $enzyme ), $tagLength ) );
                if( $bQuality == 1 ) { # get quality scores
                    my @arr = @{$$pScores[$x]};
                    my @sec = @arr[0..($tagLength-1+length($enzyme))];
                    push( @scores, \@sec );
                }
            }
        }

        my $right_tag = substr( $ditag, length( $ditag ) - $tagLength - length( $enzyme ), $tagLength );
        $right_tag = reverse( $right_tag );
        $right_tag =~ tr/ACGTacgt/TGCAtgca/;

        my $ignoreTag2 = 0;
        foreach my $ignoreTag ( @{$this->{'protocol'}->{'IGNORE_TAGS'}} ) {
            if( $ignoreTag eq $right_tag ) { $ignoreTag2 = 1; next; }
        }
        if( $ignoreTag2 == 0 ) {
            if( $right_tag =~ /^[ACGT]+$/ ) { # must contain only A,C,G,T
                push( @tags, $right_tag );
                if( $bQuality == 1 ) { # get quality scores
                    my @arr = @{$$pScores[$x]};
                    my @sec = @arr[(length($ditag)-$tagLength-length($enzyme))..(length($ditag)-1)];
                    @sec = reverse( @sec );
                    push( @scores, \@sec );
                }
            }
        }
    }

    if( $bInternal && $this->is_caching_enabled() ) {
        push( @{$this->{'tags'}}, @tags );
    }

    if( $bQuality == 1 ) {
        return ( \@tags, \@scores );
    }

    return \@tags;

}

#######################################################
sub extract_tagcounts {
#######################################################
=pod

=head2 extract_tagcounts [\@ditags], [\@scores]

Extracts valid tags from ditags and returns a hashref
containing tag sequences paired with their respective
counts.

B<Arguments>

I<\@ditags> (optional)

  An arrayref containing ditag sequences. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used.

I<\@scores> (optional)

  An array containing score arrays corresponding
  to the specified ditag sequences.

B<Returns>

  A hashref where the tag sequence is paried
  with its observed number.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @reads = ( 'ACGTAGACATAGACAAGAGATATAGA',
                'GATAGACAAAGGAAGATTACAAGATTAT' );

  my @ditags = @{$sage->extract_ditags( \@reads )};

  # get tag counts
  my %counts = %{$sage->extract_tagcounts( \@ditags )};

  # print tag counts
  map { print join( "\t", $_, $counts{$_} ) . "\n" } keys %counts; 

=cut

    my $this = shift;
    my $pDitags = shift || $this->extract_ditags();
    my $pScores = shift;
    
    my %counts;
 
    my ( $pTags, $pTagScores ) = $this->extract_tags( $pDitags, $pScores );
    map { $counts{$_}++ } @$pTags;
    
    return \%counts;

}

#######################################################
sub get_ditag_base_distribution {
#######################################################
=pod

=head2 get_ditag_base_distribution [\@ditags], [$minLength], [$maxLength]

Calculates the distribution of bases at each position and both
orientations of a set of ditags.  This distribution is used
for calculating the 'expected' nucleotide count when determining
extra bases using get_extra_base_calculation.

For example:

  CATGAAACCGTATGTAGAGAGGGACACATG
  CATGTAGACAGATAGACACAGATACCATG

has a distribution of:

          +---------------+---------------+
          |    forward    |    reverse    |
    +-----+---+---+---+---+---+---+---+---+
    | pos | A | C | G | T | A | C | G | T |
    +-----+---+---+---+---+---+---+---+---+
    |  0  | 0 | 2 | 0 | 0 | 0 | 2 | 0 | 0 |
    |  1  | 2 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
    |  2  | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 2 |
    |  3  | 0 | 0 | 2 | 0 | 0 | 0 | 2 | 0 |
    |  4  | 1 | 0 | 0 | 1 | 0 | 0 | 1 | 1 |
    |  5  | 2 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
    |  6  | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 1 |
    |                    ...              |
    | 28  | 0 | 0 | 1 | 1 | 0 | 0 | 1 | 1 |
    | 29  | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 |
    +-----+---+---+---+---+---+---+---+---+

B<Arguments>

I<\@ditags> (optional)

  An arrayref containing ditag sequences. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used.

I<$minLength> (optional)

  Ignore ditags that are smaller than this minimum length.
  If the argument is not supplied, then the minimum
  ditag length for the currently selected protocol
  is used.

I<$maxLength> (optional)

  Ignore ditags that are larger than this maximum length.
  If the argument is not supplied, then the maximum
  ditag length for the currently selected protocol
  is used.

B<Returns>

  A hashref where the key is the zero-based base position 
  index, and the value is a hashref where the key is the 
  nucleotide and the value is a hashref where the key is
  either 'fwd' or 'rev' and the value is the count of that 
  nucleotide (whew!).

  i.e. $HASH{$idx}->{'A'}->{'fwd'} = 23;

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @ditags = ( 'CATGAAACCGTATGTAGAGAGGGACACATG',
                 'CATGTAGACAGATAGACACAGATACCATG' );

  my %DIST = %{$sage->get_ditag_base_distribution( \@ditags )};

  # print distribution table
  foreach my $idx ( sort { $a <=> $b } keys %DIST ) {
      print $idx . " ";
      print join( " ", map { defined( $DIST{$idx}->{$_} ) ? $DIST{$idx}->{$_} : 0 } ( 'A','C','G','T' ) );
      print "\n";
  }

=cut

    my $this = shift;
    my $pDitags = shift;
    my $minLength = shift || $this->{'protocol'}->{'MINIMUM_DITAG_LENGTH'};
    my $maxLength = shift || $this->{'protocol'}->{'MAXIMUM_DITAG_LENGTH'};

    my $bInternal = 0;
    if( !defined( $pDitags ) ) {
        $pDitags = $this->extract_ditags();
        $bInternal = 1;
    }

    if( $this->is_caching_enabled() &&
        $bInternal == 1 &&
        defined( $this->{'cache'}->{'ditag_base_distribution'}->{$minLength.'-'.$maxLength} ) ) {
        if( $DEBUG ) {
            print STDERR "Using cached ditag base distribution.\n";
        }
        return $this->{'cache'}->{'ditag_base_distribution'}->{$minLength.'-'.$maxLength};
    }

    if( $DEBUG ) {
        print STDERR "Calculating ditag base distribution...\n";
    }
        
    my %DISTRIBUTION;

    # get forward distribution
    foreach my $ditag ( @$pDitags ) {
        my @bps = split( //, $ditag );
        if( $minLength <= length( $ditag ) && length( $ditag ) <= $maxLength ) {
            for( my $i = 0; $i < length( $ditag ); $i++ ) {
                $DISTRIBUTION{$i}->{$bps[$i]}->{'fwd'}++;
            }
        }
    }

    # get reverse distribution
    foreach my $ditag ( @$pDitags ) {
        $ditag = reverse( $ditag );
        $ditag =~ tr/ACGT/TGCA/;
        my @bps = split( //, $ditag );
        if( $minLength <= length( $ditag ) && length( $ditag ) <= $maxLength ) {
            for( my $i = 0; $i < length( $ditag ); $i++ ) {
                $DISTRIBUTION{$i}->{$bps[$i]}->{'rev'}++;
            }
        }
    }
    
    if( $DEBUG ) {
        print STDERR "  Complete.\n";
    }

    if( $this->is_caching_enabled() && $bInternal == 1 ) {
        print STDERR "Caching ditag base distribution.\n";
        $this->{'cache'}->{'ditag_base_distribution'}->{$minLength.'-'.$maxLength} = \%DISTRIBUTION;
    }
    
    return \%DISTRIBUTION;
    
}

#######################################################
sub get_ditag_length_distribution {
#######################################################
=pod

=head2 get_ditag_length_distribution [\@ditags]

Calculates the distribution of ditag lengths for a set
of ditags.

For example:

  CATGAAACCGTATGTAGAGAGGGACACATG
  CATGTAGACAGATAGACACAGATACCATG
  CATGTATCGCGGCATTACGATCTAGAACATG
  CATGACGACTATATCGATAGCTAACCATG  

has a distribution of:

    +-----+---+
    | len | N |
    +-----+---+
    |  29 | 2 |
    |  30 | 1 |
    |  31 | 1 |
    +-----+---+

B<Arguments>

I<\@ditags> (optional)

  An arrayref containing ditag sequences. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used.

B<Returns>

  A hashref where the key is the ditag length, and 
  the value is the number of ditags that have this 
  length.

  i.e. $HASH{$len} = 1024;

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @ditags = ( 'CATGAAACCGTATGTAGAGAGGGACACATG',
                 'CATGTAGACAGATAGACACAGATACCATG',
                 'CATGTATCGCGGCATTACGATCTAGAACATG',
                 'CATGACGACTATATCGATAGCTAACCATG' );

  my %DIST = %{$sage->get_ditag_length_distribution( \@ditags )};

  # print distribution table
  foreach my $idx ( sort { $a <=> $b } keys %DIST ) {
      print join( "\t", $idx, $DIST{$idx} ) . "\n";
  }

=cut
 
    my $this = shift;
    my $pDitags = shift;

    my $bInternal = 0;
    if( !defined( $pDitags ) ) {
      $pDitags = $this->extract_ditags();
      $bInternal = 1;
    }

    if( $this->is_caching_enabled() &&
        $bInternal == 1 &&
        defined( $this->{'cache'}->{'distribution'} ) ) {
      return $this->{'cache'}->{'distribution'};
    }
    
    my %DIST;
    foreach my $ditag ( @$pDitags ) {
        $DIST{length($ditag)}++;
    }

    if( $this->is_caching_enabled() && $bInternal == 1 ) {
      $this->{'cache'}->{'distribution'} = \%DIST;
    }

    return \%DIST;
    
}

#######################################################
sub get_extra_base_calculation {
#######################################################
=pod

=head2 get_extra_base_calculation $tag, [\@ditags], [\%distribution]

Searches ditags for the given tag sequence, and
using these ditags (NOTE: if tag sequence is
on right side of ditag, the entire ditag is
reverse complemented so that the search tag
is always aligned along the left side)
determines if extra bases can be assigned based
on a statistically significant increase in the
incidence of a nucleotide compared to the overall
distribution of A,C,G,T at the extra position.

The p-value for each nucleotide is calculated using
a Chi-Square test. This is essentially a 2x2
contingency table containing the observed and expected
count for a nucleotide and the sum of each observed
and expected count for the remaining nucleotides.

  For example:

  +------------------+
  | nucleotide 15    |
  +-----+----+-------+
  |     | A  | C+G+T |
  +-----+----+-------+
  | exp | 55 |  144  |
  +-----+----+-------+
  | obs | 67 |  132  |
  +-----+----+-------+
  chi-square: 3.618
  p-value: 0.05715

NOTE: This method will cause an error if the current
protocol's MAXIMUM_DITAG_LENGTH parameter does not
accomodate extra nucleotides.

SDZ: Please comment on the statistics used here.  
Is Chi-Square the best choice?  The p-values
tend to get very small, very quickly - the result
is significant nucleotides at every position
(which doesn't make sense - the terminal nucleotides
belonging to the sibling tag should be randomly
distributed).  Expected behavior seems to be seen
when a p-value threshold of 0.0001 is used - but
this seems arbitrary and unnatural.  Help!

B<Arguments>

I<$tag>

  The sequence of the tag to find extra nucleotides for.
  The sequence can be 14-mer (includes anchoring enzyme
  site) or 10-mer.

I<\@ditags> (optional)

  An arrayref containing ditag sequences. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used.

I<\%distribution> (optional)

  A hashref containing the distribution of observed bases 
  for each position of the supplied set of ditags.  If 
  this isn't specified, it is created during execution.  
  A properly formatted hash can also be obtained using 
  get_ditag_base_distribution.

B<Returns>

  A hashref where the key is the nucleotide position, and 
  the value is a hashref where the keys are the nucleotide, 
  and the value is a hashref where the keys are 'exp',
  'obs', 'exp_total', 'obs_total', 'chi_square' and 
  'pvalue'.

  i.e. $HASH{$pos}->{'A'}->{'exp'} = 34;

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my $handle;
  open( $handle, "reads.fasta" );
  $sage->load_fasta_sequence( $handle );
  close( $handle );
  my %RESULT = %{$sage->get_extra_base_calculation( "AAAACTACGA" )};
  foreach my $pos ( sort { $a <=> $b} keys %RESULT ) {
    print join( " ", 
                map { sprintf("%3i",$RESULT{$pos}->{$_}->{'pvalue'})} ('A','C','G','T') } ) . "\n";
  }

=cut

    my $this = shift @_;
    my $tag = shift @_;

    # make sure this method is applicable to the current
    # protocol
    if( $this->{'protocol'}->{'MAXIMUM_DITAG_LENGTH'} <= 2 * ( length($this->get_enzyme())+$this->{'protocol'}->{'TAG_LENGTH'} ) ) {
        die( $PACKAGE . "::get_extra_base_calculation not supported for protocols where maximum ditag size does not accomodate extra nucleotides." );
    }

    my $enzyme = $this->{'enzyme'};
    my $tagLength = $this->{'protocol'}->{'TAG_LENGTH'};

    if( $tagLength + length( $enzyme ) == length( $tag ) ) {
      # user included enzyme site; cut it off
      $tag = substr( $tag, length( $enzyme ), $tagLength );
    }
    
    my $pDitags = shift @_ || $this->extract_ditags();
    my $pDistribution = shift @_ || $this->get_ditag_base_distribution( $pDitags, 2*(length($enzyme)+$tagLength)+1 );
    
    $tag = uc( $tag );
    
    my $rtag = reverse( $tag );
    $rtag =~ tr/ACGT/TGCA/;
    
    my @ditags;
    
    # CHECK FOR x-x CONCATS? EASY MARK...
    
    # HAVE TO GET POSITION BASE BIAS DISTRIBUTION

    my $fwd = 0; # number of forward ditags
    my $rev = 0; # number of reverse ditags
    
    foreach my $ditag ( @$pDitags ) {
      if( length( $ditag ) >= 2*(length($enzyme)+$tagLength)+1 ) {
        if( $ditag =~ /^$enzyme$tag/ ) {
          push( @ditags, $ditag );
          $fwd++;
          next;
        }
        if( $ditag =~ /$rtag$enzyme$/ ) {
          $ditag = reverse( $ditag );
          $ditag =~ tr/ACGT/TGCA/;
          push( @ditags, $ditag );
          $rev++;
          next;
        }
      }
    }
    
    if( $DEBUG >= 2 ) {
        foreach my $ditag ( @ditags ) { print $ditag . "\n"; }
    }
    if( $DEBUG >= 1 ) {  
      print STDERR "Found " . scalar( @ditags ) . " / " . scalar( @$pDitags ) . " ditags.\n";
    }

    my %OBS; # holds the observed counts for each base position
    foreach my $ditag ( @ditags ) {
      my @bps = split( //, $ditag );
      for( my $i = 0; $i < scalar( @bps ); $i++ ) {
        $OBS{$i}->{$bps[$i]}++;
      }
    }

    my %RESULTS;
    my %DIST = %{$pDistribution};

    foreach my $pos ( sort { $a <=> $b } keys %OBS ) {

      my $obsA = $OBS{$pos}->{'A'} || 0;
      my $obsC = $OBS{$pos}->{'C'} || 0;
      my $obsG = $OBS{$pos}->{'G'} || 0;
      my $obsT = $OBS{$pos}->{'T'} || 0;
      my $obsN = $obsA+$obsC+$obsG+$obsT;

      # calculate contributions from fwd and rev for this position
      my $expA = ( $DIST{$pos}->{'A'}->{'fwd'} || 0 * ( $fwd / ($fwd+$rev) ) ) +
                 ( $DIST{$pos}->{'A'}->{'rev'} || 0 * ( $rev / ($fwd+$rev) ) );
      my $expC = ( $DIST{$pos}->{'C'}->{'fwd'} || 0 * ( $fwd / ($fwd+$rev) ) ) +
                 ( $DIST{$pos}->{'C'}->{'rev'} || 0 * ( $rev / ($fwd+$rev) ) );
      my $expG = ( $DIST{$pos}->{'G'}->{'fwd'} || 0 * ( $fwd / ($fwd+$rev) ) ) +
                 ( $DIST{$pos}->{'G'}->{'rev'} || 0 * ( $rev / ($fwd+$rev) ) );
      my $expT = ( $DIST{$pos}->{'T'}->{'fwd'} || 0 * ( $fwd / ($fwd+$rev) ) ) +
                 ( $DIST{$pos}->{'T'}->{'rev'} || 0 * ( $rev / ($fwd+$rev) ) );
      my $expN = $expA+$expC+$expG+$expT;

      # normalize expected to same size as observed
      $expA *= $obsN / $expN;
      $expC *= $obsN / $expN;
      $expG *= $obsN / $expN;
      $expT *= $obsN / $expN;
      $expN = $expA+$expC+$expG+$expT;

      # do chi2
      my $chiA = $this->__calculateChiSquare( $obsA, $expA, $obsN-$obsA, $expN-$expA );
      my $chiC = $this->__calculateChiSquare( $obsC, $expC, $obsN-$obsC, $expN-$expC );
      my $chiG = $this->__calculateChiSquare( $obsG, $expG, $obsN-$obsG, $expN-$expG );
      my $chiT = $this->__calculateChiSquare( $obsT, $expT, $obsN-$obsT, $expN-$expT );

      # get p-value
      my $pA = Statistics::Distributions::chisqrprob( 1, $chiA );
      my $pC = Statistics::Distributions::chisqrprob( 1, $chiC );
      my $pG = Statistics::Distributions::chisqrprob( 1, $chiG );
      my $pT = Statistics::Distributions::chisqrprob( 1, $chiT );

      # store results
      $RESULTS{$pos}->{'A'}->{'obs'} = $obsA;
      $RESULTS{$pos}->{'C'}->{'obs'} = $obsC;
      $RESULTS{$pos}->{'G'}->{'obs'} = $obsG;
      $RESULTS{$pos}->{'T'}->{'obs'} = $obsT;

      $RESULTS{$pos}->{'A'}->{'exp'} = $expA;
      $RESULTS{$pos}->{'C'}->{'exp'} = $expC;
      $RESULTS{$pos}->{'G'}->{'exp'} = $expG;
      $RESULTS{$pos}->{'T'}->{'exp'} = $expT;

      $RESULTS{$pos}->{'A'}->{'obs_total'} = $obsN-$obsA;
      $RESULTS{$pos}->{'C'}->{'obs_total'} = $obsN-$obsC;
      $RESULTS{$pos}->{'G'}->{'obs_total'} = $obsN-$obsG;
      $RESULTS{$pos}->{'T'}->{'obs_total'} = $obsN-$obsT;

      $RESULTS{$pos}->{'A'}->{'exp_total'} = $expN-$expA;
      $RESULTS{$pos}->{'C'}->{'exp_total'} = $expN-$expC;
      $RESULTS{$pos}->{'G'}->{'exp_total'} = $expN-$expG;
      $RESULTS{$pos}->{'T'}->{'exp_total'} = $expN-$expT;

      $RESULTS{$pos}->{'A'}->{'chi_square'} = $chiA;
      $RESULTS{$pos}->{'C'}->{'chi_square'} = $chiC;
      $RESULTS{$pos}->{'G'}->{'chi_square'} = $chiG;
      $RESULTS{$pos}->{'T'}->{'chi_square'} = $chiT;

      $RESULTS{$pos}->{'A'}->{'pvalue'} = $pA;
      $RESULTS{$pos}->{'C'}->{'pvalue'} = $pC;
      $RESULTS{$pos}->{'G'}->{'pvalue'} = $pG;
      $RESULTS{$pos}->{'T'}->{'pvalue'} = $pT;

    }

    return \%RESULTS;
    
}

#######################################################
sub print_ditag_base_distribution {
#######################################################
=pod

=head2 print_ditag_base_distribution [\%results]

Prints a formatted report of a ditag base distribution
result hash.

An example report looks for a list of 48967 ditags looks like:

       +-------------------------------+-------------------------------+
       |            forward            |            reverse            |
  +----+-------+-------+-------+-------+-------+-------+-------+-------+
  |pos |   A   |   C   |   G   |   T   |   A   |   C   |   G   |   T   |
  +----+-------+-------+-------+-------+-------+-------+-------+-------+
  |  0 |     0 | 48967 |     0 |     0 |     0 | 48967 |     0 |     0 |
  |  1 | 48967 |     0 |     0 |     0 | 48967 |     0 |     0 |     0 |
  |  2 |     0 |     0 |     0 | 48967 |     0 |     0 |     0 | 48967 |
  |  3 |     0 |     0 | 48967 |     0 |     0 |     0 | 48967 |     0 |
  |  4 | 12060 |  8521 | 12749 | 15592 | 12129 |  8428 | 12691 | 15687 |
  |  5 | 17677 |  8356 |  8630 | 14180 | 17278 |  8540 |  8865 | 14234 |
  |  6 | 12977 | 12649 | 11366 | 11917 | 12677 | 12449 | 11766 | 12019 |
  |  7 | 12711 | 14552 | 11060 | 10602 | 12640 | 13993 | 11384 | 10798 |
  |  8 | 12323 | 11708 |  7013 | 17883 | 12266 | 11472 |  7236 | 17905 |
  |  9 | 13814 |  8741 | 13950 | 12413 | 13692 |  8844 | 14074 | 12309 |
  | 10 | 13066 | 13347 | 10770 | 11746 | 13174 | 13293 | 10748 | 11706 |
  | 11 | 17286 | 11070 |  8177 | 12386 | 17342 | 10819 |  8420 | 12339 |
  | 12 | 13285 | 10158 | 13169 | 12309 | 13130 | 10349 | 13548 | 11885 |
  | 13 | 15542 |  9996 | 11879 | 11496 | 15367 |  9924 | 12096 | 11529 |
  | 14 | 19448 |  9491 |  8369 | 11612 | 19347 |  9545 |  8438 | 11597 |
  | 15 | 11873 | 10191 | 10786 | 16078 | 11483 | 10278 | 10617 | 16550 |
  | 16 | 12270 |  9774 |  9435 | 17444 | 12384 |  9418 |  9662 | 17446 |
  | 17 | 11148 | 12906 | 10387 | 14472 | 11433 | 12841 | 10094 | 14548 |
  | 18 | 12516 | 11297 | 10297 | 14805 | 12685 | 10988 | 10343 | 14901 |
  | 19 | 11587 |  9651 | 12003 | 15675 | 11657 |  9566 | 12174 | 15533 |
  | 20 | 12249 | 12773 | 10667 | 13231 | 12312 | 12648 | 10708 | 13255 |
  | 21 | 14983 | 10189 | 10384 | 13339 | 15077 | 10089 | 10353 | 13408 |
  | 22 | 14369 |  8953 | 13275 | 12251 | 14242 |  8571 | 13761 | 12351 |
  | 23 | 11692 | 11586 | 12960 | 12630 | 11479 | 11455 | 13282 | 12698 |
  | 24 | 12807 | 10826 | 10423 | 14863 | 12727 | 10362 | 10380 | 15422 |
  | 25 | 14435 | 12502 |  8097 | 13897 | 14413 | 12437 |  8082 | 13945 |
  | 26 |  9256 | 31392 |  3499 |  4806 |  9196 | 31383 |  3560 |  4806 |
  | 27 | 25844 | 20620 |     0 |  2503 | 25844 | 20620 |     0 |  2503 |
  | 28 | 20620 |     0 |  2503 | 25844 | 20620 |     0 |  2503 | 25844 |
  | 29 |     0 |     0 | 25844 | 20620 |     0 |     0 | 25844 | 20620 |
  | 30 |     0 |     0 | 20620 |     0 |     0 |     0 | 20620 |     0 |
  +----+-------+-------+-------+-------+-------+-------+-------+-------+

B<Arguments>

I<\%results> (optional)

  A hashref of the structure supplied by
  get_ditag_base_distribution. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used to calculate the
  distribution.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @ditags = ( 'CATGAAACCGTATGTAGAGAGGGACACATG',
                 'CATGTAGACAGATAGACACAGATACCATG' );

  $sage->print_ditag_base_distribution( $sage->get_ditag_base_distribution( \@ditags ) );

=cut

  my $this = shift;
  my $pHash = shift || $this->get_ditag_base_distribution();

  my %DIST = %{$pHash};

  print "     +-------------------------------+-------------------------------+\n";
  print "     |            forward            |            reverse            |\n";
  print "+----+-------+-------+-------+-------+-------+-------+-------+-------+\n";
  print "|pos |   A   |   C   |   G   |   T   |   A   |   C   |   G   |   T   |\n";
  print "+----+-------+-------+-------+-------+-------+-------+-------+-------+\n";

  foreach my $pos ( sort { $a <=> $b } keys %DIST ) {
    print "|";
    print sprintf( "%3i", $pos );
    print " | ";
    print join( " | ", map { sprintf( "%5i", $DIST{$pos}->{$_}->{'fwd'} || 0 ) } ('A','C','G','T') );
    print " | ";
    print join( " | ", map { sprintf( "%5i", $DIST{$pos}->{$_}->{'rev'} || 0 ) } ('A','C','G','T') );
    print " |\n";
  }

  print "+----+-------+-------+-------+-------+-------+-------+-------+-------+\n";

}

#######################################################
sub print_ditag_length_distribution {
#######################################################
=pod

=head2 print_ditag_length_distribution [\%results]

Prints a formatted report of a ditag length distribution
result hash.

An example report for a list of 49205 ditags looks like:

  +------+--------+
  |  len |  count |
  +------+--------+
  |   28 |    238 |
  |   29 |   2503 |
  |   30 |  25844 |
  |   31 |  20620 |
  +------+--------+
  |total |  49205 |
  +------+--------+

B<Arguments>

I<\%results> (optional)

  A hashref of the structure supplied by
  get_ditag_length_distribution. If this
  argument is not supplied, ditags extracted from
  reads that have been added to this object with 
  add_sequence_read are used to calculate the
  distribution.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my @ditags = ( 'CATGAAACCGTATGTAGAGAGGGACACATG',
                 'CATGTAGACAGATAGACACAGATACCATG' );

  $sage->print_ditag_base_distribution( $sage->get_ditag_base_distribution( \@ditags ) );

=cut
    
    my $this = shift;
    my $pResults = shift || $this->get_ditag_length_distribution();
    
    my %DIST = %{$pResults};
    print "+------+--------+\n";
    print "|  len |  count |\n";
    print "+------+--------+\n";
    
    my $total = 0;
    foreach my $pos ( sort { $a <=> $b } keys %DIST ) {
      print "|" . join( " | ", sprintf( "%5i", $pos ), sprintf( "%6i", $DIST{$pos} ) ) . " |\n";
      $total += $DIST{$pos};
    }

    print "+------+--------+\n";
    print "|total | " . sprintf( "%6i", $total ) . " |\n";
    print "+------+--------+\n";
   
}


#######################################################
sub print_extra_base_calculation {
#######################################################
=pod

=head2 print_extra_base_calculation \%results, <$pval_cutoff>

Prints a formatted report of an extra base result
hash.

An example report looks like (the tag is 'CATGACTTTTTCAA'):

       +------------------------------------+------------------------------------+------------------------------------+------------------------------------+
       |                 A                  |                 C                  |                 G                  |                 T                  |
  +----+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+
  |pos |  obs |  exp | tobs | texp |   pval |  obs |  exp | tobs | texp |   pval |  obs |  exp | tobs | texp |   pval |  obs |  exp | tobs | texp |   pval |
  +----+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+
  |  0 |     0      0    794  48967 1.00000 |   794  48967    794  48967 1.00000 |     0      0    794  48967 1.00000 |     0      0    794  48967 1.00000 | 
  |  1 |   794  48967    794  48967 1.00000 |     0      0    794  48967 1.00000 |     0      0    794  48967 1.00000 |     0      0    794  48967 1.00000 | 
  |  2 |     0      0    794  48967 1.00000 |     0      0    794  48967 1.00000 |     0      0    794  48967 1.00000 |   794  48967    794  48967 1.00000 | 
  |  3 |     0      0    794  48967 1.00000 |     0      0    794  48967 1.00000 |   794  48967    794  48967 1.00000 |     0      0    794  48967 1.00000 | 
  |  4 |   794  12060    794  48923 0.00000 |     0   8521    794  48923 0.00000 |     0  12749    794  48923 0.00000 |     0  15592    794  48923 0.00000 | A
  |  5 |     0  17677    794  48843 0.00000 |   794   8356    794  48843 0.00000 |     0   8630    794  48843 0.00000 |     0  14180    794  48843 0.00000 | C
  |  6 |     0  12977    794  48909 0.00000 |     0  12649    794  48909 0.00000 |     0  11366    794  48909 0.00000 |   794  11917    794  48909 0.00000 | T
  |  7 |     0  12711    794  48925 0.00000 |     0  14552    794  48925 0.00000 |     0  11060    794  48925 0.00000 |   794  10602    794  48925 0.00000 | T
  |  8 |     0  12323    794  48927 0.00000 |     0  11708    794  48927 0.00000 |     0   7013    794  48927 0.00000 |   794  17883    794  48927 0.00000 | T
  |  9 |     0  13814    794  48919 0.00000 |     0   8741    794  48919 0.00000 |     0  13950    794  48919 0.00000 |   794  12413    794  48919 0.00000 | T
  | 10 |     0  13066    794  48929 0.00000 |     0  13347    794  48929 0.00000 |     0  10770    794  48929 0.00000 |   794  11746    794  48929 0.00000 | T
  | 11 |     0  17286    794  48919 0.00000 |   794  11070    794  48919 0.00000 |     0   8177    794  48919 0.00000 |     0  12386    794  48919 0.00000 | C
  | 12 |   794  13285    794  48921 0.00000 |     0  10158    794  48921 0.00000 |     0  13169    794  48921 0.00000 |     0  12309    794  48921 0.00000 | A
  | 13 |   794  15542    794  48913 0.00000 |     0   9996    794  48913 0.00000 |     0  11879    794  48913 0.00000 |     0  11496    794  48913 0.00000 | A
  | 14 |   768  19448    793  48921 0.00000 |     5   9491    793  48921 0.00000 |     9   8369    793  48921 0.00000 |    11  11612    793  48921 0.00000 | A
  | 15 |   284  11873    794  48929 0.00000 |   135  10191    794  48929 0.00069 |   216  10786    794  48929 0.00000 |   159  16078    794  48929 0.00000 | AG
  | 16 |   212  12270    794  48923 0.13956 |   213   9774    794  48923 0.00000 |   198   9435    794  48923 0.00000 |   171  17444    794  48923 0.00000 | CG
  | 17 |   151  11148    794  48913 0.00079 |   230  12906    794  48913 0.01589 |   178  10387    794  48913 0.29035 |   235  14472    794  48913 0.99248 | 
  | 18 |   190  12516    794  48915 0.12272 |   221  11297    794  48915 0.00001 |   172  10297    794  48915 0.59363 |   211  14805    794  48915 0.00044 | C
  | 19 |   189  11587    793  48917 0.91485 |   177   9651    793  48917 0.01775 |   202  12003    793  48917 0.38779 |   225  15675    793  48917 0.00028 | 
  | 20 |   192  12249    794  48921 0.41673 |   233  12773    794  48921 0.00213 |   166  10667    794  48921 0.40376 |   203  13231    794  48921 0.17060 | 
  | 21 |   231  14983    794  48895 0.13224 |   180  10189    794  48895 0.09725 |   191  10384    794  48895 0.00914 |   192  13339    794  48895 0.00401 | 
  | 22 |   219  14369    792  48849 0.09370 |   166   8953    792  48849 0.01544 |   219  13275    792  48849 0.65444 |   188  12251    792  48849 0.21110 | 
  | 23 |   175  11692    792  48869 0.09870 |   214  11586    792  48869 0.00222 |   178  12960    792  48869 0.00022 |   225  12630    792  48869 0.01578 | 
  | 24 |   222  12807    794  48919 0.10078 |   193  10826    794  48919 0.04471 |   146  10423    794  48919 0.00866 |   233  14863    794  48919 0.31392 | 
  | 25 |   232  14435    792  48931 0.81419 |   214  12502    792  48931 0.18065 |   138   8097    792  48931 0.42601 |   208  13897    792  48931 0.04213 | 
  | 26 |   121   9256    794  48953 0.00106 |   574  31392    794  48953 1.00000 |    40   3499    794  48953 0.01337 |    59   4806    794  48953 0.01282 | 
  | 27 |   466  25844    794  48967 1.00000 |   287  20620    794  48967 0.00000 |     0      0    794  48967 1.00000 |    41   2503    794  48967 0.94376 | 
  | 28 |   287  20620    794  48967 0.00000 |     0      0    794  48967 1.00000 |    41   2503    794  48967 0.94376 |   466  25844    794  48967 1.00000 | 
  | 29 |     0      0    753  46465 1.00000 |     0      0    753  46465 1.00000 |   466  25844    753  46465 1.00000 |   287  20620    753  46465 0.00000 | 
  | 30 |     0      0    287  20621 1.00000 |     0      0    287  20621 1.00000 |   287  20620    287  20621 1.00000 |     0      0    287  20621 1.00000 | 
  +----+------------------------------------+------------------------------------+------------------------------------+------------------------------------+

Using this report, one could reasonably conclude that the 15th base of a 14-mer SAGE tag is an 'A'.

B<Arguments>

I<\%results>

  A hashref of the structure created using
  get_extra_base_calculation.

I<$pval_cutoff> (optional)

  The report prints out statistically significant
  nucleotides at the right of the table.  This argument
  is the p-value that is used to consider if the nucleotide
  is present at a statistically significant number.

  The default value is 0.0001.

B<Usage>

  my $sage = Bio::SAGE::DataProcessing->new();
  my $handle;
  open( $handle, "reads.fasta" );
  $sage->load_fasta_sequence( $handle );
  close( $handle );
  my %results = %{$sage->get_extra_base_calculation( "AAAACTACGA" )};
  $sage->print_extra_base_calculation ( \%results );

=cut

    my $this = shift;
    my $pHash = shift;
    my $pval_cutoff = shift || 0.0001;

    my %RESULTS = %{$pHash};

    # print header
    print "     +------------------------------------+------------------------------------+------------------------------------+------------------------------------+\n";
    print "     |                 A                  |                 C                  |                 G                  |                 T                  |\n";
    print "+----+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+\n";
    print "|pos |  obs |  exp | tobs | texp |   pval |  obs |  exp | tobs | texp |   pval |  obs |  exp | tobs | texp |   pval |  obs |  exp | tobs | texp |   pval |\n";
    print "+----+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+------+------+------+------+--------+\n";

    foreach my $pos ( sort { $a <=> $b } keys %RESULTS ) {

      my @bps = ( 'A','C','G','T' );
      print "|";
      print join( " | ", sprintf("%3i",$pos),
                         map {
                           join( " ", sprintf("%5i",$RESULTS{$pos}->{$_}->{'obs'}),
                                      sprintf("%6i",$RESULTS{$pos}->{$_}->{'exp'}),
                                      sprintf("%6i",$RESULTS{$pos}->{$_}->{'obs_total'}),
                                      sprintf("%6i",$RESULTS{$pos}->{$_}->{'exp_total'}),
                                      #sprintf("%3.2f", $RESULTS{$pos}->{$_}->{'chi_square'}),
                                      sprintf("%1.5f", $RESULTS{$pos}->{$_}->{'pvalue'}) )
                         } @bps );

      print " | ";
      print join( "", map { if( $RESULTS{$pos}->{$_}->{'pvalue'} < $pval_cutoff &&
                                $RESULTS{$pos}->{$_}->{'obs'}/( $RESULTS{$pos}->{$_}->{'obs_total'} == 0 ? 1 : $RESULTS{$pos}->{$_}->{'obs_total'} ) >= 
                                $RESULTS{$pos}->{$_}->{'exp'}/( $RESULTS{$pos}->{$_}->{'exp_total'} == 0 ? 1 : $RESULTS{$pos}->{$_}->{'exp_total'} ) ) {
                              $_;
                            }
                          } @bps ) . "\n";

    }

    # print footer
    print "+----+------------------------------------+------------------------------------+------------------------------------+------------------------------------+\n";
        
}

sub __calculateChiSquare {

    my $this = shift;

    my $obs = shift; # observed count
    my $exp = shift; # expected count
    my $obsTotal = shift; # total observations
    my $expTotal = shift || 1; # total expecteds

    # normalize expected
    $exp *= $obsTotal / $expTotal;
    $expTotal = $obsTotal;

    #print join( " | ", $obs, $exp, $obsTotal, $expTotal );

    my $chi1 = $obs - $exp; $chi1 *= $chi1; $chi1 /= $exp || 1;
    my $chi2 = ( $obsTotal - $obs ) - ( $expTotal - $exp ); $chi2 *= $chi2; 
    $chi2 /= ( $expTotal == $exp ? 1 : ( $expTotal - $exp ) );

    return ( $chi1 + $chi2 );

}

sub __convertPhred2Probability {

  my $this = shift;
  my $phred = shift;

  $phred /= 10;
  $phred *= -1;

  return 10 ** $phred;

}

sub __calculateProbability {

  my $this = shift;
  my $pScores = shift;

  my $p = 1.0;
  foreach my $score ( @$pScores ) {
    $p *= ( 1.0 - $this->__convertPhred2Probability( $score ) );
  }

  return $p

}

sub __clearCache {

  my $this = shift;
  
  %{$this->{'cache'}} = ();

}

1;

__END__

=pod

=head1 COPYRIGHT

Copyright(c)2004 Scott Zuyderduyn <scottz@bccrc.ca>. All rights reserved.

This program is free software; you can redistribute it and/or modify it 
under the same terms as Perl itself.

=head1 AUTHOR

Scott Zuyderduyn <scottz@bccrc.ca>
BC Cancer Research Centre

=head1 VERSION

  1.00

=head1 SEE ALSO

Statistics::Distributions(1).

=head1 TODO

  Add more debugging messages.

=cut
