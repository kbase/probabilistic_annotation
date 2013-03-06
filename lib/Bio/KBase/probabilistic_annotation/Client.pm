package Bio::KBase::probabilistic_annotation::Client;

use JSON::RPC::Client;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

Bio::KBase::probabilistic_annotation::Client

=head1 DESCRIPTION





=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => Bio::KBase::probabilistic_annotation::Client::RpcClient->new,
	url => $url,
    };


    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 annotation_probabilities

  $return = $obj->annotation_probabilities($annotation_probabilities_input)

=over 4

=item Parameter and return types

=begin html

<pre>
$annotation_probabilities_input is an annotation_probabilities_input
$return is a ProbabilisticAnnotation
annotation_probabilities_input is a reference to a hash where the following keys are defined:
	probanno_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	genomeObj has a value which is a GenomeObject
workspace_id is a string
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
	workspace has a value which is a workspace_id
ProbAnnoFeature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	alternative_functions has a value which is a reference to a list where each element is an alt_func

</pre>

=end html

=begin text

$annotation_probabilities_input is an annotation_probabilities_input
$return is a ProbabilisticAnnotation
annotation_probabilities_input is a reference to a hash where the following keys are defined:
	probanno_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	genomeObj has a value which is a GenomeObject
workspace_id is a string
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
	workspace has a value which is a workspace_id
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
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotation_probabilities (received $n, expecting 1)");
    }
    {
	my($annotation_probabilities_input) = @args;

	my @_bad_arguments;
        (ref($annotation_probabilities_input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"annotation_probabilities_input\" (value was \"$annotation_probabilities_input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotation_probabilities:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotation_probabilities');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "ProbabilisticAnnotation.annotation_probabilities",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'annotation_probabilities',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotation_probabilities",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotation_probabilities',
				       );
    }
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
	probanno_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	genome has a value which is a genome_id
workspace_id is a string
probanno_id is a string
genome_id is a string
ProbabilisticAnnotation is a reference to a hash where the following keys are defined:
	id has a value which is a probanno_id
	genome has a value which is a genome_id
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
	workspace has a value which is a workspace_id
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
	probanno_workspace has a value which is a workspace_id
	probanno has a value which is a probanno_id
	genome has a value which is a genome_id
workspace_id is a string
probanno_id is a string
genome_id is a string
ProbabilisticAnnotation is a reference to a hash where the following keys are defined:
	id has a value which is a probanno_id
	genome has a value which is a genome_id
	featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
	workspace has a value which is a workspace_id
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
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotation_probabilities_id (received $n, expecting 1)");
    }
    {
	my($annotation_probabilities_ids_input) = @args;

	my @_bad_arguments;
        (ref($annotation_probabilities_ids_input) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"annotation_probabilities_ids_input\" (value was \"$annotation_probabilities_ids_input\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotation_probabilities_id:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotation_probabilities_id');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "ProbabilisticAnnotation.annotation_probabilities_id",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'annotation_probabilities_id',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotation_probabilities_id",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotation_probabilities_id',
				       );
    }
}



sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, {
        method => "ProbabilisticAnnotation.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'annotation_probabilities_id',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method annotation_probabilities_id",
            status_line => $self->{client}->status_line,
            method_name => 'annotation_probabilities_id',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for Bio::KBase::probabilistic_annotation::Client\n";
    }
    if ($sMajor == 0) {
        warn "Bio::KBase::probabilistic_annotation::Client version is $svr_version. API subject to change.\n";
    }
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
workspace has a value which is a workspace_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a probanno_id
genome has a value which is a genome_id
featureAlternativeFunctions has a value which is a reference to a list where each element is a ProbAnnoFeature
workspace has a value which is a workspace_id


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
probanno_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
genomeObj has a value which is a GenomeObject

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
probanno_workspace has a value which is a workspace_id
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
probanno_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
genome has a value which is a genome_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
probanno_workspace has a value which is a workspace_id
probanno has a value which is a probanno_id
genome has a value which is a genome_id


=end text

=back



=cut

package Bio::KBase::probabilistic_annotation::Client::RpcClient;
use base 'JSON::RPC::Client';

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $obj) = @_;
    my $result;

    if ($uri =~ /\?/) {
       $result = $self->_get($uri);
    }
    else {
        Carp::croak "not hashref." unless (ref $obj eq 'HASH');
        $result = $self->_post($uri, $obj);
    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        $obj->{id} = $self->id if (defined $self->id);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
