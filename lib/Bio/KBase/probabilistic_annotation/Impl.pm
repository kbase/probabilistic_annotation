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



=head2 annotate

  $jobid = $obj->annotate($input)

=over 4

=item Parameter and return types

=begin html

<pre>
$input is an annotate_params
$jobid is a string
annotate_params is a reference to a hash where the following keys are defined:
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

</pre>

=end html

=begin text

$input is an annotate_params
$jobid is a string
annotate_params is a reference to a hash where the following keys are defined:
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


=end text



=item Description

This function, rather than using an already-loaded genome object, loads a genome from the specified workspace
before running the probabilistic annotation algorithm.

=back

=cut

sub annotate
{
    my $self = shift;
    my($input) = @_;

    my @_bad_arguments;
    (ref($input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"input\" (value was \"$input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($jobid);
    #BEGIN annotate
    
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
    
    # Run the python code.
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
	
    #END annotate
    my @_bad_returns;
    (!ref($jobid)) or push(@_bad_returns, "Invalid type for return variable \"jobid\" (value was \"$jobid\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate');
    }
    return($jobid);
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
    my $cmdline = "probanno-normalize --model '".$input->{model}."' --modelws ".$input->{model_workspace};
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



=head2 genome_id

=over 4



=item Description

A string identifier for a genome.


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



=item Description

A string identifier for a feature.


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



=head2 FunctionProbability

=over 4



=item Description

Annotation probability for an alternative function

        string function - the name of the functional role being annotated to the feature
        float probability - the probability that the functional role is associated with the feature
        string functionMD5 - hash to let us know if anything has changed


=item Definition

=begin html

<pre>
a reference to a list containing 3 items:
0: (function) a string
1: (probability) a float
2: (functionMD5) a string

</pre>

=end html

=begin text

a reference to a list containing 3 items:
0: (function) a string
1: (probability) a float
2: (functionMD5) a string


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



=head2 annotate_params

=over 4



=item Description

Input parameters for the "annotate" function.

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
