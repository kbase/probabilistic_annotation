#!/usr/bin/env perl

# Need a license here

use strict;
use warnings;
use Getopt::Long::Descriptive;
use Bio::KBase::probabilistic_annotation::Client;
use Bio::KBase::probabilistic_annotation::Helpers qw(get_probanno_client);
use Bio::KBase::workspaceService::Helpers qw(auth workspace printObjectMeta);

my $manpage =
"
NAME
      pa-annotate -- generate probabilistic annotation for a genome

SYNOPSIS
      pa-annotate <Genome ID> <ProbAnno ID> [OPTIONS]

DESCRIPTION
      Generate a probabilistic annotation for a genome.
      
      Options:
      -e, --showerror    Show any errors in execution
      --genomews         ID of workspace where Genome object is stored
      -h, --help         Display this help message, ignore all arguments
      -o, --overwrite    Overwrite existing ProbAnno object with same name
      -v, --verbose      Print verbose messages
      -w, --probannows   ID of workspace where ProbAnno object is saved

EXAMPLES
      Annotate:
      > pa-annotate kb\|g.0 kb.g.0
      
AUTHORS
      Matt Benedict, Mike Mundy
";

# Define usage and options.
my $primaryArgs = [ "Genome ID", "ProbAnno ID" ];
my ( $opt, $usage ) = describe_options(
    'pa-annotate <' . join( "> <", @{$primaryArgs} ) . '> %o',
    [ 'probannows:s', 'ID of workspace where ProbAnno object is saved', { "default" => workspace() } ],
    [ 'genomews:s', 'ID of workspace where Genome object is stored', { "default" => "KBaseCDMGenomes" } ],
    [ 'overwrite:o', "Set as 1 to overwrite existing ProbAnno object with same name", { "default" => 0 } ],
    [ 'showerror|e', 'Set as 1 to show any errors in execution', { "default" => 0 } ],
    [ 'verbose|v', 'Set as 1 to print verbose messages', { "default" => 0 } ],
    [ 'help|h', 'Show help text' ],
    [ 'usage|?', 'Show usage information' ]
    );
if ( defined( $opt->{help} ) ) {
    print $manpage;
    exit 0;
}
if (defined($opt->{usage})) {
	print $usage;
	exit 0;
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
my $client = get_probanno_client();

# Define translation from options to function parameters.
my $translation = {
    "Genome ID"   => "genome",
    "ProbAnno ID" => "probanno",
    genomews      => "genome_workspace",
    probannows    => "probanno_workspace",
    overwrite     => "overwrite"
};

# Instantiate parameters for function.
my $params = { auth => auth(), };
foreach my $key ( keys( %{$translation} ) ) {
    if ( defined( $opt->{$key} ) ) {
		$params->{ $translation->{$key} } = $opt->{$key};
    }
}
print "overwrite=".$opt->{overwrite}."\n";

# Call the function.
my $output = $client->annotation_probabilities_id($params);
if (!defined($output)) {
	print "Probabilistic annotation failed!\n";
	exit 1;
} else {
	print "Probabilistic annotation successfully generated in workspace:\n";
	printObjectMeta($output);
}
exit 0;
