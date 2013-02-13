#!/usr/bin/perl
#
# This function automatically looks for a JSON file for the genome
# (this has to be gotten from IRIS for now - package of KBase for linux
# doesn't have the necessary tools yet!) and runs the probability calculation
# service with the genome object as input.
#
use strict;
use Bio::KBase::probabilistic_annotation::Impl
#use JSON;
#use File::Spec;

# Load up a test genome JSON file
# I only do this to get a genome object - in teh future this won't be necessary
# because I'll either be able to call annotate_genome to get a genome object
# or to get a genome object for a particular ID in the CDM.
my $genomeid = $ARGV[0] or die "Usage: test_impl.pl [genomeID]";
#my $jsonFileName = File::Spec->catfile("${genomeid}", "${genomeid}.json");
#open(FILE, "<${jsonFileName}") or die "Unable to open expected input JSON file ${jsonFileName}";
#my $jsonString = join("", <FILE>);
#close(FILE);
#my $genomeObject = decode_json $jsonString;

# Create a new object that will be computing the probabilities for me.
my $annoteObject = ProbabilisticAnnotation->new();

# Call the function - this function will first write the genome object to a JSON file (in our case,
# the same JSON file as was used to make the genome object), then call the pyhton script that uses it
# to do its calculations and modifies it by adding probabilities, and finally will read back and return
# an object containing the probabilities (genomeTO)
my $genomeTO = $annoteObject->annotation_probabilities_id($genomeid);

print "Done!\n";
