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
      pa-gendata -- generate data files for probabilistic annotation

SYNOPSIS
      pa-gendata <Folder path> [OPTIONS]

DESCRIPTION
      Generate data files for probabilistic annotation of genomes.
      
      Options:
      -d, --deleteonly   Delete data files but do not regenerate
      -h, --help         Display this help message, ignore all arguments
      -r, --regenerate   Regenerate data files even if they exist
      -v, --verbose      Print verbose messages
      
      Note that --deleteonly and --regenerate are mutually exclusive.

EXAMPLES
      Annotate:
      > pa-gendata data
      
AUTHORS
      Matt Benedict, Mike Mundy
";

# Define usage and options.
my $primaryArgs = [ "Folder path" ];
my ( $opt, $usage ) = describe_options(
    'pa-gendata <' . join( "> <", @{$primaryArgs} ) . '> %o',
    [ 'deleteonly|d', 'Set as 1 to delete but not regenerate data files', { "default" => 0 } ],
    [ 'regenerate|r', 'Set as 1 to delete and regenerate data files', { "default" => 0 } ],
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
    "Folder path" => "folder",
    deleteonly    => "delete_only",
    regenerate    => "regenerate",
    verbose       => "verbose"
};

# Instantiate parameters for function.
my $params = { auth => auth(), };
foreach my $key ( keys( %{$translation} ) ) {
    if ( defined( $opt->{$key} ) ) {
		$params->{ $translation->{$key} } = $opt->{$key};
    }
}

# Validate parameters.
if ($params->{delete_only} && $params->{regenerate}) {
	print "You cannot specify both --regenerate and --deleteonly\n";
	print $manpage;
	exit(1);
}

# Call the function.
my $success = $client->generate_data($params);
if (!defined($success) || $success == 0) {
	print "Generating data failed!\n";
	exit 1;
} else {
	print "Data files successfully generated in specified folder\n";
}
exit 0;
