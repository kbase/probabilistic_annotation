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
      pa-normalize -- generate compound weights for a probability model

SYNOPSIS
      pa-normalize <Model ID> [OPTIONS]

DESCRIPTION
      Generate compound (metabolite) weights for a probability model.
      
      Options:
      -d, --debug        Keep intermediate files for debug purposes
      -e, --showerror    Show any errors in execution
      -h, --help         Display this help message, ignore all arguments
      -v, --verbose      Print verbose messages
      -w, --modelws      ID of workspace where Model object is saved

EXAMPLES
      Annotate:
      > pa-normalize kb\|g.8622
      
AUTHORS
      Matt Benedict, Mike Mundy
";

# Define usage and options.
my $primaryArgs = [ "Model ID" ];
my ( $opt, $usage ) = describe_options(
    'pa-normalize <' . join( "> <", @{$primaryArgs} ) . '> %o',
    [ 'modelws|w:s', 'ID of workspace where Model object is saved', { "default" => workspace() } ],
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
    "Model ID"    => "model",
    modelws       => "model_workspace",
    debug         => "debug"
};

# Instantiate parameters for function.
my $params = { auth => auth(), };
foreach my $key ( keys( %{$translation} ) ) {
    if ( defined( $opt->{$key} ) ) {
		$params->{ $translation->{$key} } = $opt->{$key};
    }
}

# Call the function.
my $output = $client->normalize($params);
if (!defined($output)) {
	print "Normalization of model failed!\n";
	exit 1;
} else {
	print "Probability model successfully normalized\n";
}
exit 0;
