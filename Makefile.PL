# *%) $Id: Makefile.PL,v 1.1.1.1 2004/05/24 02:52:55 scottz Exp $

use ExtUtils::MakeMaker;

#$VERSION = "1.00";

WriteMakefile(
   'DISTNAME'  => "Bio-SAGE-DataProcessing",
   'NAME'      => "Bio::SAGE::DataProcessing",
   'AUTHOR'    => 'Scott Zuyderduyn <scottz@bccrc.ca>',
   'ABSTRACT'  => 'Module to process raw data generated from SAGE libraries.',
   'VERSION_FROM' => 'lib/Bio/SAGE/DataProcessing.pm',
#   'VERSION'   => $VERSION,
   'PM'        => { 'lib/Bio/SAGE/DataProcessing.pm' => '$(INST_LIBDIR)/DataProcessing.pm', },
   'PREREQ_PM' => { 'Statistics::Distributions' => 1.02, },
);

