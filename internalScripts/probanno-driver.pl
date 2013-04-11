#!/usr/bin/perl
#
# This function automatically looks for a JSON file for the genome
# (this has to be gotten from IRIS for now - package of KBase for linux
# doesn't have the necessary tools yet!) and runs the probability calculation
# service with the genome object as input.
#
use strict;
use Bio::KBase::AuthToken;
use Bio::KBase::probabilistic_annotation::Impl;
use Bio::KBase::workspaceService::Impl;
use Bio::KBase::workspaceService::Helpers qw(auth get_ws_client);

#my $token = Bio::KBase::AuthToken->new(token => auth());
#print $token->user_id();

# Load up a test genome JSON file
# I only do this to get a genome object - in teh future this won't be necessary
# because I'll either be able to call annotate_genome to get a genome object
# or to get a genome object for a particular ID in the CDM.
my $usage = "Usage: test_impl.pl [genomeID] [genome_workspace_id] [probannoID] [probanno_workspace_id]";
my $genomeid = $ARGV[0] or die $usage;
my $genome_workspace_id = $ARGV[1] or die $usage;
my $probanno_id = $ARGV[2] or die $usage;
my $probanno_workspace_id = $ARGV[3] or die $usage;

# Check workspace permissions for currently logged-in user
my $serv = get_ws_client();
my $genome_metadata = $serv->get_workspacemeta( { auth => auth(), 
						  workspace => $genome_workspace_id } );
my $probanno_metadata = $serv->get_workspacemeta( { auth => auth(), 
						    workspace => $probanno_workspace_id } );
#print ${$probanno_metadata}[4];
if ( ${$genome_metadata}[4] eq "n" ) {
    die "At least read (r) permission is required for genome workspace";
}
if ( ${$probanno_metadata}[4] eq "n" || ${$probanno_metadata}[4] eq "r" ) {
    die "Either w or a permission is required for probanno workspace";
}

my $annotation_probabilities_id_input = { "genome_workspace"   => $genome_workspace_id,
					  "probanno"           => $probanno_id,
					  "genome"             => $genomeid };

# Create a new object that will be computing the probabilities for me.
my $annoteObject = Bio::KBase::probabilistic_annotation::Impl->new();

# Call the function - this function will first write the genome object to a JSON file (in our case,
# the same JSON file as was used to make the genome object), then call the pyhton script that uses it
# to do its calculations and modifies it by adding probabilities, and finally will read back and return
# an object containing the probabilities (genomeTO)
my $probanno_object = $annoteObject->annotation_probabilities_id($annotation_probabilities_id_input);

# Save the results to the desired probanno_workspace
my $servercommand = "save_object";

my $params = { "id" => $probanno_id,
	       "workspace" => $probanno_workspace_id,
	       "auth" => auth(),
	       "type" => "ProbAnno",
	       "data" => $probanno_object
};

my $savedObject = $serv->$servercommand($params);

print "Done!\n";
