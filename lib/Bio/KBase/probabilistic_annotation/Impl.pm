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
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternative_functions has a value which is a reference to a list where each element is an alt_func

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
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternative_functions has a value which is a reference to a list where each element is an alt_func


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
    my $workspacefolder = "/kb/deployment/data/probabilistic_annotation";
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

  $return = $obj->annotation_probabilities_id($annotation_probabilities_ids_input)

=over 4

=item Parameter and return types

=begin html

<pre>
$annotation_probabilities_ids_input is an annotation_probabilities_ids_input
$return is a ProbabilisticAnnotation
annotation_probabilities_ids_input is a reference to a hash where the following keys are defined:
	genome_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	genome has a value which is a genome_id
workspace_id is a string
probanno_id is a string
genome_id is a string
ProbabilisticAnnotation is a reference to a hash where the following keys are defined:
	id has a value which is a probanno_id
	genome has a value which is a genome_id
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternative_functions has a value which is a reference to a list where each element is an alt_func
feature_id is a string
alt_func is a reference to a list containing 2 items:
	0: (function) a string
	1: (probability) a float

</pre>

=end html

=begin text

$annotation_probabilities_ids_input is an annotation_probabilities_ids_input
$return is a ProbabilisticAnnotation
annotation_probabilities_ids_input is a reference to a hash where the following keys are defined:
	genome_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	genome has a value which is a genome_id
workspace_id is a string
probanno_id is a string
genome_id is a string
ProbabilisticAnnotation is a reference to a hash where the following keys are defined:
	id has a value which is a probanno_id
	genome has a value which is a genome_id
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternative_functions has a value which is a reference to a list where each element is an alt_func
feature_id is a string
alt_func is a reference to a list containing 2 items:
	0: (function) a string
	1: (probability) a float


=end text



=item Description

This function, rather than using an already-loaded genome object, loads a genome from the specified workspace
before running the probabilistic annotation algorithm.

=back

=cut

sub annotation_probabilities_id
{
    my $self = shift;
    my($annotation_probabilities_ids_input) = @_;

    my @_bad_arguments;
    (ref($annotation_probabilities_ids_input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"annotation_probabilities_ids_input\" (value was \"$annotation_probabilities_ids_input\")");
    my $badref = ref($annotation_probabilities_ids_input);
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotation_probabilities_id - got $badref :\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities_id');
    }

    my $ctx = $Bio::KBase::probabilistic_annotation::Server::CallContext;
    my($return);
    #BEGIN annotation_probabilities_id

    # I'm not entirely sure if this is the right approach but attempting to use get_object
    # directly from the workspace object failed - it told me it couldn't connect to the server.
    my $serv = get_ws_client();
    my $servercommand = "get_object";

    my $workspacefolder = "/kb/deployment/data/probabilistic_annotation";
    mkdir($workspacefolder);

    # Note - we relegate importing the new object to a workspace to another script (particularly the wrapper script for this API call). 
    # This script just returns the probanno object to be imported.
    my $genome_id = $annotation_probabilities_ids_input->{"genome"};
    my $genome_workspace_id = $annotation_probabilities_ids_input->{"genome_workspace"};
    my $probanno_id = $annotation_probabilities_ids_input->{"probanno"};

    # The output stuff should also go to a standard place.
    # Both of these should be changed to point at the workspace's data folder if there is one...
    my $outputdir = File::Spec->catdir("$workspacefolder", "$genome_id");
    mkdir($outputdir);

    my $jsonFileName = File::Spec->catfile("$outputdir", "${genome_id}.json");
   
    # Pull genome from the specified workspace
    my $objParams = { "id" => $genome_id,
		      "workspace" => $genome_workspace_id,
		      "type" => "Genome",
		      auth => auth() };

    # WHY does the workspace unwrap the data when calling kbws but not when calling it through this???
    # I'm pretty sure I must be doing this wrong... as always.
    my $genomeObject = $serv->$servercommand($objParams)->{"data"};
    
    my $annotation_probabilities_input = { "probanno"           => $probanno_id,
					   "genomeObj"          => $genomeObject };

    $return = $self->annotation_probabilities($annotation_probabilities_input);

    #END annotation_probabilities_id
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotation_probabilities_id:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities_id');
    }
    return($return);
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



=head2 annotationProbability

=over 4



=item Description

Data structures to hold a single annotation probability for a single gene

feature_id feature - feature the annotation is associated with
string function - the name of the functional role being annotated to the feature
float probability - the probability that the functional role is associated with the feature


=item Definition

=begin html

<pre>
a reference to a list containing 3 items:
0: (feature) a feature_id
1: (function) a string
2: (probability) a float

</pre>

=end html

=begin text

a reference to a list containing 3 items:
0: (feature) a feature_id
1: (function) a string
2: (probability) a float


=end text

=back



=head2 probanno_id

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



=head2 ProbAnnoFeature

=over 4



=item Description

Object to carry alternative functions for each feature
    
feature_id id
ID of the feature. Required.
    
string function
Primary annotated function of the feature in the genome annotation. Required.
    
list<alt_func> alternative_functions
List of tuples containing alternative functions and probabilities. Required.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a feature_id
alternative_functions has a value which is a reference to a list where each element is an alt_func

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a feature_id
alternative_functions has a value which is a reference to a list where each element is an alt_func


=end text

=back



=head2 workspace_id

=over 4



=item Description

Object to carry alternative functions and probabilities for genes in a genome    

        probanno_id id - ID of the probabilistic annotation object. Required.    
        genome_id genome - ID of the genome the probabilistic annotation was built for. Required.
        list<ProbAnnoFeature> featureAlternativeFunctions - List of ProbAnnoFeature objects holding alternative functions for features. Required.    
        workspace - ID of the workspace from which the genome ID was taken. Required.


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



=head2 ProbabilisticAnnotation

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature


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



=head2 annotation_probabilities_ids_input

=over 4



=item Description

I will implement getting the genome from the workspace later. For now its a headache that is too much
beyond all the other stuff I have to change here.

This genome_id is a CDM genome_id.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
genome has a value which is a genome_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
genome has a value which is a genome_id


=end text

=back



=cut

1;
