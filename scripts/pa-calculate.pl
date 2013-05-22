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
      pa-calculate <ProbAnno ID> <RxnProb ID> [OPTIONS]

DESCRIPTION
      Generate reaction probabilities from a probabilistic annotation.
      
      Options:
      -d, --debug           Keep intermediate files for debug purposes
      -e, --showerror       Show any errors in execution
      -h, --help            Display this help message, ignore all arguments
      -o, --overwrite       Overwrite existing Model object with same name
      --probannows          ID of workspace where ProbAnno object is stored
      -t, --templatemodel   Template model ID
      -m, --templatemodelws Template model workspace
      -v, --verbose         Print verbose messages
      -w, --rxnprobws        Workspace for reaction probability object

EXAMPLES
      Annotate:
      > pa-calculate kb\|g.0 kb\|g.0
      
AUTHORS
      Matt Benedict, Mike Mundy
";

# Define usage and options.
my $primaryArgs = [ "ProbAnno ID", "RxnProb ID" ];
my ( $opt, $usage ) = describe_options(
    'pa-calculate <' . join( "> <", @{$primaryArgs} ) . '> %o',
    [ 'probannows|w=s', 'ID of workspace where ProbAnno object is stored', { "default" => workspace() } ],
    [ 'rxnprobws|w=s', 'ID of workspace to save the output', { 'default' => workspace() } ],
    [ 'templatemodel|t=s', "template model", { "default" => undef } ],
    [ 'templatemodelws|m=s', "template model workspace", { "default" => undef } ],
    [ 'debug|d:i', "Keep intermediate files for debug purposes", { "default" => 0 } ],
    [ 'showerror|e:i', 'Show any errors in execution', { "default" => 0 } ],
    [ 'verbose|v:i', 'Print verbose messages', { "default" => 0 } ],
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
		print STDERR "Required arguments are missing\n".$usage;
		exit 1;
    }
}

# Create a client object.
my $client = get_probanno_client();

# Define translation from options to function parameters.
my $translation = {
    "ProbAnno ID" => "probanno",
    "RxnProb ID"  => "outputid",
    "rxnprobws"    => "outputws",
    probannows    => "probanno_workspace",
    debug         => "debug",
    verbose       => "verbose",
    templatemodel => "template_model",
    templatemodelws => "template_model_workspace"
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
