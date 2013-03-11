#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long::Descriptive;
use Bio::KBase::probabilistic_annotation::Client;
use Bio::KBase::workspaceService::Helpers qw(auth workspace workspaceURL);

# Need to make this a helper function.
#my $probannoURL = "http://localhost:7073/";
my $probannoURL = "http://140.221.84.212:7073";

# Define usage and options.
my $primaryArgs = [ "Genome ID", "ProbAnno ID" ];
my ( $opt, $usage ) = describe_options(
    'pa-annotate <' . join( "> <", @{$primaryArgs} ) . '> %o',
    [ 'workspace|w:s', 'Workspace to load probabilistic annotation into', { "default" => workspace() } ],
    [ 'genomews:s', 'Workspace with genome', { "default" => "KBaseCDMGenomes" } ],
    [ 'overwrite:o', "Overwrite existing probablistic annotation with same name", { "default" => 0 } ],
    [ 'showerror|e', 'Set as 1 to show any errors in execution', { "default" => 0 } ],
    [ 'verbose|v', 'Print verbose messages', { "default" => 0 } ],
    [ 'help|h|?',  'Print this usage information' ]
    );
if ( defined( $opt->{help} ) ) {
    print $usage;
    exit;
}

# Process primary arguments.
foreach my $arg ( @{$primaryArgs} ) {
    $opt->{$arg} = shift @ARGV;
    if ( !defined( $opt->{$arg} ) ) {
		print $usage;
		exit;
    }
}

# Create a client object.
my $client = Bio::KBase::probabilistic_annotation::Client->new($probannoURL);

# Define translation from options to annotate() parameters.
my $translation = {
    "Genome ID" => "genome",
    "ProbAnno ID" => "probanno",
    genomews => "genome_workspace",
    workspace => "workspace",
    overwrite => "overwrite"
};

# Instantiate parameters for annotate() function.
my $params = { auth => auth(), };
foreach my $key ( keys( %{$translation} ) ) {
    if ( defined( $opt->{$key} ) ) {
		$params->{ $translation->{$key} } = $opt->{$key};
    }
}

print ref($params);
print "\n";

# Call the function.
my $output = $client->annotation_probabilities_id($params);
if (!defined($output)) {
	print "Probabilistic annotation failed!\n";
} else {
	print "Probabilistic annotation succeeded!\n";
#	print $output;
}
