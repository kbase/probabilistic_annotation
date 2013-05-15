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
      pa-calculate -- generate probability model for a genome

SYNOPSIS
      pa-calculate <ProbAnno ID> <Model ID> [OPTIONS]

DESCRIPTION
      Generate a probability model for a genome.
      
      Options:
      -d, --debug        Keep intermediate files for debug purposes
      -e, --showerror    Show any errors in execution
      -h, --help         Display this help message, ignore all arguments
      -o, --overwrite    Overwrite existing Model object with same name
      --probannows       ID of workspace where ProbAnno object is stored
      -v, --verbose      Print verbose messages
      -w, --modelws      ID of workspace where Model object is saved

EXAMPLES
      Annotate:
      > pa-calculate kb\|g.0 kb\|g.0
      
AUTHORS
      Matt Benedict, Mike Mundy
";

# Define usage and options.
my $primaryArgs = [ "ProbAnno ID" ];
my ( $opt, $usage ) = describe_options(
    'pa-calculate <' . join( "> <", @{$primaryArgs} ) . '> %o',
    [ 'probannows:s', 'ID of workspace where ProbAnno object is stored', { "default" => workspace() } ],
    [ 'overwrite|o', "Set as 1 to overwrite existing Model object with same name", { "default" => 0 } ],
    [ 'debug|d', "Set as 1 to keep intermediate files for debug purposes", { "default" => 0 } ],
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
    "ProbAnno ID" => "probanno",
    probannows    => "probanno_workspace",
    debug         => "debug",
    verbose       => "verbose"
};

# Instantiate parameters for function.
my $params = { auth => auth(), };
foreach my $key ( keys( %{$translation} ) ) {
    if ( defined( $opt->{$key} ) ) {
		$params->{ $translation->{$key} } = $opt->{$key};
    }
}

# Call the function.
my $output = $client->calculate($params);
if (!defined($output)) {
	print "Calculating reactions failed!\n";
	exit 1;
} else {
	foreach my $rxnprob ( @{$output}) {
		print join("\t", $rxnprob->[0], $rxnprob->[1], $rxnprob->[2], $rxnprob->[3], $rxnprob->[4])."\n";
	}
}
exit 0;
