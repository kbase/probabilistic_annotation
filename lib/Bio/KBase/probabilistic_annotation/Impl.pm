package Bio::KBase::probabilistic_annotation::Impl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

ProbabilisticAnnotation

=head1 DESCRIPTION



=cut

#BEGIN_HEADER
use JSON::XS;
use File::Spec;
use Bio::KBase::workspaceService::Impl;
use Bio::KBase::workspaceService::Helpers qw(auth get_ws_client);
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 annotation_probabilities

  $return = $obj->annotation_probabilities($annotation_probabilities_input)

=over 4

=item Parameter and return types

=begin html

<pre>
$annotation_probabilities_input is an annotation_probabilities_input
$return is a ProbabilisticAnnotation
annotation_probabilities_input is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	genomeObj has a value which is a GenomeObject
probanno_id is a string
GenomeObject is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	contigs has a value which is a reference to a list where each element is a contig
	features has a value which is a reference to a list where each element is a feature
genome_id is a string
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
contig_id is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	alternative_functions has a value which is a reference to a list where each element is an alt_func
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	annotations has a value which is a reference to a list where each element is an annotation
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
alt_func is a reference to a list containing 2 items:
	0: (function) a string
	1: (probability) a float
annotation is a reference to a list containing 3 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
ProbabilisticAnnotation is a reference to a hash where the following keys are defined:
	id has a value which is a probanno_id
	genome has a value which is a genome_id
	genome_workspace has a value which is a workspace_id
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
	rolesetProbabilities has a value which is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a FunctionProbability
	skippedFeatures has a value which is a reference to a list where each element is a feature_id
workspace_id is a string
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternativeFunctions has a value which is a reference to a list where each element is a FunctionProbability
FunctionProbability is a reference to a list containing 2 items:
	0: (function) a string
	1: (probability) a float

</pre>

=end html

=begin text

$annotation_probabilities_input is an annotation_probabilities_input
$return is a ProbabilisticAnnotation
annotation_probabilities_input is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	genomeObj has a value which is a GenomeObject
probanno_id is a string
GenomeObject is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	contigs has a value which is a reference to a list where each element is a contig
	features has a value which is a reference to a list where each element is a feature
genome_id is a string
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
contig_id is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	alternative_functions has a value which is a reference to a list where each element is an alt_func
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	annotations has a value which is a reference to a list where each element is an annotation
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
alt_func is a reference to a list containing 2 items:
	0: (function) a string
	1: (probability) a float
annotation is a reference to a list containing 3 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
ProbabilisticAnnotation is a reference to a hash where the following keys are defined:
	id has a value which is a probanno_id
	genome has a value which is a genome_id
	genome_workspace has a value which is a workspace_id
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
	rolesetProbabilities has a value which is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a FunctionProbability
	skippedFeatures has a value which is a reference to a list where each element is a feature_id
workspace_id is a string
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternativeFunctions has a value which is a reference to a list where each element is a FunctionProbability
FunctionProbability is a reference to a list containing 2 items:
	0: (function) a string
	1: (probability) a float


=end text



=item Description

Given a genome object populated with genes and annotations, this function adds
potential alternative functions with probabilities

=back

=cut

sub annotation_probabilities
{
    my $self = shift;
    my($annotation_probabilities_input) = @_;

    my @_bad_arguments;
    (ref($annotation_probabilities_input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"annotation_probabilities_input\" (value was \"$annotation_probabilities_input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotation_probabilities:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($return);
    #BEGIN annotation_probabilities

    # Get the genome ID from the genome object and using it to create a home for our output files.
    my $genomeObject = $annotation_probabilities_input->{"genomeObj"};
    my $probanno_id = $annotation_probabilities_input->{"probanno"};
    my $genome_id = $genomeObject->{"id"};

#    $genome_id =~ s/\|/_/;

    # I've hard-coded this for now because I can't get the references to work through the symlink that the
    # dev_container generates.
    #
    # This is on my TODO list to fix.
#    my $workspacefolder = "/kb/deployment/data/probabilistic_annotation";
    my $workspacefolder = "data";
    mkdir($workspacefolder);

    # The output stuff should also go to a standard place.
    # Both of these should be changed to point at the workspace's data folder if there is one...
    my $outputdir = File::Spec->catdir("$workspacefolder", "$genome_id");
    mkdir($outputdir);

    # Make the JSON string and dump it to a file.
    # ASCII - I'm not sure what the best way to deal with this is but I don't want it to just die
    # if a non-UTF8 character is encountered in the roles (which has happened to me from time to time...)
    my $JSON_STRING = JSON::XS->new->ascii->pretty->encode($genomeObject);
    my $outfile = File::Spec->catfile("$outputdir","${genome_id}.json");

    open(FILE, ">$outfile") or die "Unable to create file ${outfile} to which to dump the provided genome object to annotation_probabilities";
    print FILE $JSON_STRING;
    close(FILE);

    # System call to probability calculator
    # This script must be in the PATH or we will fail.
    # It is a wrapped-up version of Probability_calculation_frontend.py
    my $status = system("Probability_calculation_frontend", "-f", "$workspacefolder", "$genome_id", "$probanno_id");
    if ( ($status >>= 8) != 0 ) {
	die "Probability calculator failed.";
    }

    # Read the new JSON file (for the new probability object type)
    my $infile = File::Spec->catfile("${outputdir}","${genome_id}_prob.json");
    print "${infile}\n";
    open(FILE, "<$infile") or die "Unable to read file ${infile} which should contain the probabilities from annotation_probabilities";
    $JSON_STRING = join("", <FILE>); 
    close(FILE);
    $return = JSON::XS->new->ascii->decode($JSON_STRING);
    
    #END annotation_probabilities
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotation_probabilities:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities');
    }
    return($return);
}




=head2 annotation_probabilities_id

  $output = $obj->annotation_probabilities_id($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is an annotation_probabilities_ids_params
$output is an object_metadata
annotation_probabilities_ids_params is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	genome_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	overwrite has a value which is a bool
	debug has a value which is a bool
	auth has a value which is a string
genome_id is a string
workspace_id is a string
probanno_id is a string
bool is an int
object_metadata is a reference to a list containing 11 items:
	0: (id) an object_id
	1: (type) an object_type
	2: (moddate) a timestamp
	3: (instance) an int
	4: (command) a string
	5: (lastmodifier) a username
	6: (owner) a username
	7: (workspace) a workspace_id
	8: (ref) a workspace_ref
	9: (chsum) a string
	10: (metadata) a reference to a hash where the key is a string and the value is a string
object_id is a string
object_type is a string
timestamp is a string
username is a string
workspace_ref is a string

</pre>

=end html

=begin text

$input is an annotation_probabilities_ids_params
$output is an object_metadata
annotation_probabilities_ids_params is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	genome_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	overwrite has a value which is a bool
	debug has a value which is a bool
	auth has a value which is a string
genome_id is a string
workspace_id is a string
probanno_id is a string
bool is an int
object_metadata is a reference to a list containing 11 items:
	0: (id) an object_id
	1: (type) an object_type
	2: (moddate) a timestamp
	3: (instance) an int
	4: (command) a string
	5: (lastmodifier) a username
	6: (owner) a username
	7: (workspace) a workspace_id
	8: (ref) a workspace_ref
	9: (chsum) a string
	10: (metadata) a reference to a hash where the key is a string and the value is a string
object_id is a string
object_type is a string
timestamp is a string
username is a string
workspace_ref is a string


=end text



=item Description

This function, rather than using an already-loaded genome object, loads a genome from the specified workspace
before running the probabilistic annotation algorithm.

=back

=cut

sub annotation_probabilities_id
{
    my $self = shift;
    my($input) = @_;

    my @_bad_arguments;
    (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"input\" (value was \"$input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotation_probabilities_id:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities_id');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($output);
    #BEGIN annotation_probabilities_id

    # Build command line to bridge to Python script.
    my $cmdline = "probanno-annotate --genome '".$input->{genome}."' --genomews ".$input->{genome_workspace};
    $cmdline .= " --probanno '".$input->{probanno}."' --probannows ".$input->{probanno_workspace};
    $cmdline .= " --auth '".$input->{auth}."'";
    if ($input->{overwrite}) {
    	$cmdline .= " --overwrite 1";
    }
    if ($input->{debug}) {
    	$cmdline .= " --debug 1";
    }
    if ($input->{verbose}) {
    	$cmdline .= " --verbose 1";
    }
    print STDERR $cmdline."\n";
    my $status = system($cmdline);
    if ($status == -1) {
     	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-annotate failed to execute",
    		method_name => 'annotation_probabilities_id');	   	
    }
    elsif ($status & 127) {
     	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-annotate failed with signal ".$status&127,
    		method_name => 'annotation_probabilities_id');	   	    	
    }
    elsif ($status >> 8 == 0) {
    	# This is lame but needed until I switch to pure python server.
	    my $ws_client = get_ws_client();
	    my $getobjmeta_params = {
	       id         => $input->{probanno},
	       workspace  => $input->{probanno_workspace},
	       type       => "ProbAnno",
	       auth       => $input->{auth}
       };
	   $output = $ws_client->get_objectmeta($getobjmeta_params);
    }
    else {
    	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-annotate failed with return code ".$status>>8,
    		method_name => 'annotation_probabilities_id');	
    }
	
    #END annotation_probabilities_id
    my @_bad_returns;
    (ref($output) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotation_probabilities_id:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities_id');
    }
    return($output);
}




=head2 calculate

  $output = $obj->calculate($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a calculate_params
$output is an object_metadata
calculate_params is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	model has a value which is a model_id
	model_workspace has a value which is a workspace_id
	overwrite has a value which is a bool
	debug has a value which is a bool
	auth has a value which is a string
probanno_id is a string
workspace_id is a string
model_id is a string
bool is an int
object_metadata is a reference to a list containing 11 items:
	0: (id) an object_id
	1: (type) an object_type
	2: (moddate) a timestamp
	3: (instance) an int
	4: (command) a string
	5: (lastmodifier) a username
	6: (owner) a username
	7: (workspace) a workspace_id
	8: (ref) a workspace_ref
	9: (chsum) a string
	10: (metadata) a reference to a hash where the key is a string and the value is a string
object_id is a string
object_type is a string
timestamp is a string
username is a string
workspace_ref is a string

</pre>

=end html

=begin text

$input is a calculate_params
$output is an object_metadata
calculate_params is a reference to a hash where the following keys are defined:
	probanno has a value which is a probanno_id
	probanno_workspace has a value which is a workspace_id
	model has a value which is a model_id
	model_workspace has a value which is a workspace_id
	overwrite has a value which is a bool
	debug has a value which is a bool
	auth has a value which is a string
probanno_id is a string
workspace_id is a string
model_id is a string
bool is an int
object_metadata is a reference to a list containing 11 items:
	0: (id) an object_id
	1: (type) an object_type
	2: (moddate) a timestamp
	3: (instance) an int
	4: (command) a string
	5: (lastmodifier) a username
	6: (owner) a username
	7: (workspace) a workspace_id
	8: (ref) a workspace_ref
	9: (chsum) a string
	10: (metadata) a reference to a hash where the key is a string and the value is a string
object_id is a string
object_type is a string
timestamp is a string
username is a string
workspace_ref is a string


=end text



=item Description



=back

=cut

sub calculate
{
    my $self = shift;
    my($input) = @_;

    my @_bad_arguments;
    (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"input\" (value was \"$input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to calculate:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'calculate');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($output);
    #BEGIN calculate

    # Build command line to bridge to Python script.
    my $cmdline = "probanno-calculate --probanno '".$input->{probanno}."' --probannows ".$input->{probanno_workspace};
    $cmdline .= " --model '".$input->{model}."' --modelws ".$input->{model_workspace};
    $cmdline .= " --auth '".$input->{auth}."'";
    if ($input->{overwrite}) {
    	$cmdline .= " --overwrite 1";
    }
    if ($input->{debug}) {
    	$cmdline .= " --debug 1";
    }
    if ($input->{verbose}) {
    	$cmdline .= " --verbose 1";
    }

	# Run the python code.
    my $status = system($cmdline);
    if ($status == -1) {
     	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-calculate failed to execute",
    		method_name => 'calculate');	   	
    }
    elsif ($status & 127) {
     	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-calculate failed with signal ".$status&127,
    		method_name => 'calculate');	   	    	
    }
    elsif ($status >> 8 == 0) {
    	# This is lame but needed until I switch to pure python server.
	    my $ws_client = get_ws_client();
	    my $getobjmeta_params = {
	       id         => $input->{model},
	       workspace  => $input->{model_workspace},
	       type       => "Model",
	       auth       => $input->{auth}
       };
	   $output = $ws_client->get_objectmeta($getobjmeta_params);
    }
    else {
    	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-calculate failed with return code ".$status>>8,
    		method_name => 'calculate');	
    }

    #END calculate
    my @_bad_returns;
    (ref($output) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to calculate:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'calculate');
    }
    return($output);
}




=head2 normalize

  $success = $obj->normalize($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a normalize_params
$success is a bool
normalize_params is a reference to a hash where the following keys are defined:
	model has a value which is a model_id
	model_workspace has a value which is a workspace_id
	debug has a value which is a bool
	auth has a value which is a string
model_id is a string
workspace_id is a string
bool is an int

</pre>

=end html

=begin text

$input is a normalize_params
$success is a bool
normalize_params is a reference to a hash where the following keys are defined:
	model has a value which is a model_id
	model_workspace has a value which is a workspace_id
	debug has a value which is a bool
	auth has a value which is a string
model_id is a string
workspace_id is a string
bool is an int


=end text



=item Description



=back

=cut

sub normalize
{
    my $self = shift;
    my($input) = @_;

    my @_bad_arguments;
    (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"input\" (value was \"$input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to normalize:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'normalize');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($success);
    #BEGIN normalize
    
    # Build command line to bridge to Python script.
    my $cmdline = "probanno-normalize --model ".$input->{model}." -modelws ".$input->{model_workspace};
    $cmdline .= " --auth '".$input->{auth}."'";
    if ($input->{debug}) {
    	$cmdline .= " --debug 1";
    }
    
    # Run the python code.
    my $status = system($cmdline);
    if ($status == -1) {
     	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-normalize failed to execute",
    		method_name => 'normalize');	   	
    }
    elsif ($status & 127) {
     	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-normalize failed with signal ".$status&127,
    		method_name => 'normalize');	   	    	
    }
    elsif ($status >> 8 == 0) {
    	$success = 1;
    }
    else {
    	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => "probanno-normalize failed with return code ".$status>>8,
    		method_name => 'normalize');	
    }
    #END normalize
    my @_bad_returns;
    (!ref($success)) or push(@_bad_returns, "Invalid type for return variable \"success\" (value was \"$success\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to normalize:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'normalize');
    }
    return($success);
}




=head2 generate_data

  $success = $obj->generate_data($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is a generate_data_params
$success is a bool
generate_data_params is a reference to a hash where the following keys are defined:
	folder has a value which is a string
	regenerate has a value which is a bool
	delete_only has a value which is a bool
	verbose has a value which is a bool
bool is an int

</pre>

=end html

=begin text

$input is a generate_data_params
$success is a bool
generate_data_params is a reference to a hash where the following keys are defined:
	folder has a value which is a string
	regenerate has a value which is a bool
	delete_only has a value which is a bool
	verbose has a value which is a bool
bool is an int


=end text



=item Description



=back

=cut

sub generate_data
{
    my $self = shift;
    my($input) = @_;

    my @_bad_arguments;
    (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"input\" (value was \"$input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to generate_data:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'generate_data');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($success);
    #BEGIN generate_data
    
    $success = 1;
    
    # Build command line to bridge to Python script.
    my $cmdline = "probanno-gendata -f ".$input->{folder};
    if ($input->{delete_only}) {
    	$cmdline .= " -d";
    }
    if ($input->{regenerate}) {
    	$cmdline .= " -r";
    }
    if ($input->{verbose}) {
    	$cmdline .= " -v";
    }
    my $status = system($cmdline);
    if ( ($status >>= 8) != 0 ) {
		$success = 0;
    }
    
    #END generate_data
    my @_bad_returns;
    (!ref($success)) or push(@_bad_returns, "Invalid type for return variable \"success\" (value was \"$success\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to generate_data:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'generate_data');
    }
    return($success);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 bool

=over 4



=item Description

Indicates true or false values (false <= 0, true >=1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 probanno_id

=over 4



=item Description

A string identifier for a probabilistic annotation object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_id

=over 4



=item Description

A string identifier for a workspace. Any string consisting of alphanumeric characters and "-" is acceptable.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 object_type

=over 4



=item Description

A string indicating the type of an object stored in a workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 model_id

=over 4



=item Description

A string identifier for a model object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 object_id

=over 4



=item Description

ID of an object stored in the workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 username

=over 4



=item Description

Login name of KBase user account to which permissions for workspaces are mapped


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 timestamp

=over 4



=item Description

Exact time for workspace operations. e.g. 2012-12-17T23:24:06


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_ref

=over 4



=item Description

A permanent reference to an object in a workspace.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 object_metadata

=over 4



=item Description

Meta data associated with an object stored in a workspace.

        object_id id - ID of the object assigned by the user or retreived from the IDserver (e.g. kb|g.0)
        object_type type - type of the object (e.g. Genome)
        timestamp moddate - date when the object was modified by the user (e.g. 2012-12-17T23:24:06)
        int instance - instance of the object, which is equal to the number of times the user has overwritten the object
        timestamp date_created - time at which the alignment was built/loaded in seconds since the epoch
        string command - name of the command last used to modify or create the object
        username lastmodifier - name of the user who last modified the object
        username owner - name of the user who owns (who created) this object
        workspace_id workspace - ID of the workspace in which the object is currently stored
        workspace_ref ref - a 36 character ID that provides permanent undeniable access to this specific instance of this object
        string chsum - checksum of the associated data object
        mapping<string,string> metadata - custom metadata entered for data object during save operation


=item Definition

=begin html

<pre>
a reference to a list containing 11 items:
0: (id) an object_id
1: (type) an object_type
2: (moddate) a timestamp
3: (instance) an int
4: (command) a string
5: (lastmodifier) a username
6: (owner) a username
7: (workspace) a workspace_id
8: (ref) a workspace_ref
9: (chsum) a string
10: (metadata) a reference to a hash where the key is a string and the value is a string

</pre>

=end html

=begin text

a reference to a list containing 11 items:
0: (id) an object_id
1: (type) an object_type
2: (moddate) a timestamp
3: (instance) an int
4: (command) a string
5: (lastmodifier) a username
6: (owner) a username
7: (workspace) a workspace_id
8: (ref) a workspace_ref
9: (chsum) a string
10: (metadata) a reference to a hash where the key is a string and the value is a string


=end text

=back



=head2 md5

=over 4



=item Description

********************************************************************************
    Genome object (from CDM) type definition
    ********************************************************************************


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 md5s

=over 4



=item Definition

=begin html

<pre>
a reference to a list where each element is a md5
</pre>

=end html

=begin text

a reference to a list where each element is a md5

=end text

=back



=head2 genome_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 contig_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 feature_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 feature_type

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 region_of_dna

=over 4



=item Description

A region of DNA is maintained as a tuple of four components:

        the contig
        the beginning position (from 1)
        the strand
        the length

        We often speak of "a region".  By "location", we mean a sequence
        of regions from the same genome (perhaps from distinct contigs).


=item Definition

=begin html

<pre>
a reference to a list containing 4 items:
0: a contig_id
1: (begin) an int
2: (strand) a string
3: (length) an int

</pre>

=end html

=begin text

a reference to a list containing 4 items:
0: a contig_id
1: (begin) an int
2: (strand) a string
3: (length) an int


=end text

=back



=head2 location

=over 4



=item Description

a "location" refers to a sequence of regions


=item Definition

=begin html

<pre>
a reference to a list where each element is a region_of_dna
</pre>

=end html

=begin text

a reference to a list where each element is a region_of_dna

=end text

=back



=head2 annotation

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 3 items:
0: (comment) a string
1: (annotator) a string
2: (annotation_time) an int

</pre>

=end html

=begin text

a reference to a list containing 3 items:
0: (comment) a string
1: (annotator) a string
2: (annotation_time) an int


=end text

=back



=head2 alt_func

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: (function) a string
1: (probability) a float

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: (function) a string
1: (probability) a float


=end text

=back



=head2 feature

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a feature_id
location has a value which is a location
type has a value which is a feature_type
function has a value which is a string
alternative_functions has a value which is a reference to a list where each element is an alt_func
protein_translation has a value which is a string
aliases has a value which is a reference to a list where each element is a string
annotations has a value which is a reference to a list where each element is an annotation

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a feature_id
location has a value which is a location
type has a value which is a feature_type
function has a value which is a string
alternative_functions has a value which is a reference to a list where each element is an alt_func
protein_translation has a value which is a string
aliases has a value which is a reference to a list where each element is a string
annotations has a value which is a reference to a list where each element is an annotation


=end text

=back



=head2 contig

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a contig_id
dna has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a contig_id
dna has a value which is a string


=end text

=back



=head2 GenomeObject

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a genome_id
scientific_name has a value which is a string
domain has a value which is a string
genetic_code has a value which is an int
source has a value which is a string
source_id has a value which is a string
contigs has a value which is a reference to a list where each element is a contig
features has a value which is a reference to a list where each element is a feature

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a genome_id
scientific_name has a value which is a string
domain has a value which is a string
genetic_code has a value which is an int
source has a value which is a string
source_id has a value which is a string
contigs has a value which is a reference to a list where each element is a contig
features has a value which is a reference to a list where each element is a feature


=end text

=back



=head2 FunctionProbability

=over 4



=item Description

Annotation probability for an alternative function

        string function - the name of the functional role being annotated to the feature
        float probability - the probability that the functional role is associated with the feature


=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: (function) a string
1: (probability) a float

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: (function) a string
1: (probability) a float


=end text

=back



=head2 ProbAnnoFeature

=over 4



=item Description

Alternative functions for each feature

    feature_id id - ID of feature the annotation is associated with 
    list<FunctionProbability> alternativeFunctions - list of alternative functions and probabilities


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a feature_id
alternativeFunctions has a value which is a reference to a list where each element is a FunctionProbability

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a feature_id
alternativeFunctions has a value which is a reference to a list where each element is a FunctionProbability


=end text

=back



=head2 ProbabilisticAnnotation

=over 4



=item Description

Object to carry alternative functions and probabilities for genes in a genome    

        probanno_id id - ID of the probabilistic annotation object    
        genome_id genome - ID of the genome the probabilistic annotation was built for
        workspace_id genome_workspace - ID of the workspace containing genome
        list<ProbAnnoFeature> featureAlternativeFunctions - list of ProbAnnoFeature objects holding alternative functions for features
        mapping<feature_id feature, list<FunctionProbability>> rolesetProbabilities - mapping of features to list of FunctionProbability objects
        list<feature_id> skippedFeatures - list of features in genome with no probability


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
rolesetProbabilities has a value which is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a FunctionProbability
skippedFeatures has a value which is a reference to a list where each element is a feature_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
rolesetProbabilities has a value which is a reference to a hash where the key is a feature_id and the value is a reference to a list where each element is a FunctionProbability
skippedFeatures has a value which is a reference to a list where each element is a feature_id


=end text

=back



=head2 annotation_probabilities_input

=over 4



=item Description

********************************************** 
                     Function definitions
    ***********************************************


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
genomeObj has a value which is a GenomeObject

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
genomeObj has a value which is a GenomeObject


=end text

=back



=head2 annotation_probabilities_ids_params

=over 4



=item Description

Input parameters for the "annotation_probabilities_ids" function.

       genome_id genome - ID of Genome object
       workspace_id genome_workspace - ID of workspace where Genome object is stored
       probanno_id probanno - ID of ProbAnno object
       workspace_id probanno_workspace - ID workspace where ProbAnno object is saved
       bool overwrite - True to overwrite existing ProbAnno object with same name
       bool debug - True to keep intermediate files for debug purposes
       string auth - Authentication token of KBase user


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
overwrite has a value which is a bool
debug has a value which is a bool
auth has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome has a value which is a genome_id
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
overwrite has a value which is a bool
debug has a value which is a bool
auth has a value which is a string


=end text

=back



=head2 calculate_params

=over 4



=item Description

Input parameters for the "calculate" function.

            probanno_id probanno - ID of ProbAnno object
            workspace_id probanno_workspace - ID of workspace where ProbAnno object is stored
            model_id model - ID of Model object
            workspace_id model_workspace - ID of workspace where Model object is saved   
            bool overwrite - True to overwrite existing ProbAnno object with same name
            bool debug - True to keep intermediate files for debug purposes
            string auth - Authentication token of KBase user


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
model has a value which is a model_id
model_workspace has a value which is a workspace_id
overwrite has a value which is a bool
debug has a value which is a bool
auth has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
probanno has a value which is a probanno_id
probanno_workspace has a value which is a workspace_id
model has a value which is a model_id
model_workspace has a value which is a workspace_id
overwrite has a value which is a bool
debug has a value which is a bool
auth has a value which is a string


=end text

=back



=head2 normalize_params

=over 4



=item Description

Input parameters for the "normalize" function.

            model_id model - ID of Model object
            workspace_id model_workspace - ID of workspace where Model object is saved   
            bool debug - True to keep intermediate files for debug purposes
            string auth - Authentication token of KBase user


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
model has a value which is a model_id
model_workspace has a value which is a workspace_id
debug has a value which is a bool
auth has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
model has a value which is a model_id
model_workspace has a value which is a workspace_id
debug has a value which is a bool
auth has a value which is a string


=end text

=back



=head2 generate_data_params

=over 4



=item Description

Input parameters for the "generate_data" function.

        string folder - Path to folder for generated data files
        bool regenerate - True to delete and regenerate existing data files
        bool delete_only - True to only delete existing data files
        bool verbose - True to enable verbose output


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
folder has a value which is a string
regenerate has a value which is a bool
delete_only has a value which is a bool
verbose has a value which is a bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
folder has a value which is a string
regenerate has a value which is a bool
delete_only has a value which is a bool
verbose has a value which is a bool


=end text

=back



=cut

1;
