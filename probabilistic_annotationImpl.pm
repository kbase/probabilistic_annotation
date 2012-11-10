package probabilistic_annotationImpl;
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

  $return = $obj->annotation_probabilities($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
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
	1: an int
	2: a string
	3: an int
feature_type is a string
alt_func is a reference to a list containing 2 items:
	0: a string
	1: a float
annotation is a reference to a list containing 3 items:
	0: a string
	1: a string
	2: an int

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
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
	1: an int
	2: a string
	3: an int
feature_type is a string
alt_func is a reference to a list containing 2 items:
	0: a string
	1: a float
annotation is a reference to a list containing 3 items:
	0: a string
	1: a string
	2: an int


=end text



=item Description

Given a genome object populated with genes and annotations, this function adds
potential alternative functions with probabilities

=back

=cut

sub annotation_probabilities
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotation_probabilities:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities');
    }

    my $ctx = $probabilistic_annotationServer::CallContext;
    my($return);
    #BEGIN annotation_probabilities

    # Get the genome ID from the genome object and using it to create a home for our output files.
    my $genomeId = $genomeTO->{"id"};
    mkdir($genomeId);
    # Make the JSON string and dump it to a file.
    # ASCII - I'm not sure what the best way to deal with this is but I don't want it to just die
    # if a non-UTF8 character is encountered in the roles (which has happened to me from time to time...)
    my $JSON_STRING = JSON::XS->new->ascii->pretty->encode($genomeTO);
    my $outfile = File::Spec->catfile("${genomeId}","${genomeId}.json");
    open(FILE, ">$outfile") or die "Unable to create file ${outfile} to which to dump the provided genome object to annotation_probabilities";
    print FILE $JSON_STRING;
    close(FILE);
    # System call to probability calculator

    # CHRIS: The folder where the workspace function ultimately dumps the files should be
    # placed in this variable
    #
    # Fow now I'm using a relative path "data" just so you can see what it does.
    # It will create the folder if it doesn't exist and dump files to it.
    #
    # The expected names of each file are listed in PYTHON_GLOBALS.py. If possible please
    # have the workspace extractor save any temporary files with those names. If this is not
    # possible I'll change that.
    my $workspacefolder = "data";
    my $status = system("python", "Probability_calculation_frontend.py", "-f", "$workspacefolder", "$genomeId");
    if ( ($status >>= 8) != 0 ) {
	die "Probability calculator failed.";
    }
    # Read the new JSON file (for the new probability object type)
    my $infile = File::Spec->catfile("${genomeId}","${genomeId}_prob.json");
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

  $return = $obj->annotation_probabilities_id($genome_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_id is a genome_id
$return is a genomeTO
genome_id is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	contigs has a value which is a reference to a list where each element is a contig
	features has a value which is a reference to a list where each element is a feature
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
	1: an int
	2: a string
	3: an int
feature_type is a string
alt_func is a reference to a list containing 2 items:
	0: a string
	1: a float
annotation is a reference to a list containing 3 items:
	0: a string
	1: a string
	2: an int

</pre>

=end html

=begin text

$genome_id is a genome_id
$return is a genomeTO
genome_id is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	contigs has a value which is a reference to a list where each element is a contig
	features has a value which is a reference to a list where each element is a feature
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
	1: an int
	2: a string
	3: an int
feature_type is a string
alt_func is a reference to a list containing 2 items:
	0: a string
	1: a float
annotation is a reference to a list containing 3 items:
	0: a string
	1: a string
	2: an int


=end text



=item Description

This does the same thing except it takes a genome ID and attempts to search for
the specified ID in the central store before throwing an exception

=back

=cut

sub annotation_probabilities_id
{
    my $self = shift;
    my($genome_id) = @_;

    my @_bad_arguments;
    (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_id\" (value was \"$genome_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotation_probabilities_id:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotation_probabilities_id');
    }

    my $ctx = $probabilistic_annotationServer::CallContext;
    my($return);
    #BEGIN annotation_probabilities_id
    # Chris said he's going to wrap the cs_to_genome function in an API routine.
    # When that is ready, it will go here:
    #my $genomeTO = [package]->cs_to_genome($genome_id);
    #
    # I've hacked this together until then so we have SOMETHING.
    my $jsonFileName = File::Spec->catfile("${genome_id}", "${genome_id}.json");
    my $jsonString;
    if ( -e $jsonFileName) {

	open(FILE, "<${jsonFileName}") or die "Unable to open input JSON file ${jsonFileName} despite it existing???";
    } else {
	open (my $file, '>', 'TEMPORARY_GENOME_JSON');
	print $file `perl cs_to_genome_MODIFIED.pl "${genome_id}"`;
	open(FILE, "<TEMPORARY_GENOME_JSON") or die "Unable to open the temporary JSON file created by cs_to_genome_MODIFIED.pl";
    }
    $jsonString = join("", <FILE>);
    close(FILE);

    my $genomeObject = decode_json $jsonString;

    $return = $self->annotation_probabilities($genomeObject);

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
1: an int
2: a string
3: an int

</pre>

=end html

=begin text

a reference to a list containing 4 items:
0: a contig_id
1: an int
2: a string
3: an int


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
0: a string
1: a string
2: an int

</pre>

=end html

=begin text

a reference to a list containing 3 items:
0: a string
1: a string
2: an int


=end text

=back



=head2 alt_func

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: a string
1: a float

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: a string
1: a float


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



=head2 genomeTO

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



=cut

1;
