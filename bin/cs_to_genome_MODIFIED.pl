
use strict;
use gjoseqlib;
use Bio::KBase::IDServer::Client;
use Bio::KBase::CDMI::Client;
use JSON::XS;
use Data::Dumper;

use Getopt::Long;

my $scientific_name;
my $domain;
my $genetic_code;
my $source;
my $source_id;
my $input_file;
my $output_file;

my $rc = GetOptions('output=s'    => \$output_file,
		    );

# MODIFIED - temporary stopgap to get a genome object with IDs as actually
# present in the KBase repo...
# This should obviously become an input option instead
# but I need to move forward with this.
my $usage = "cs_to_genome_MODIFIED [--output output-file] KBase-genome-id > genome-object";

@ARGV == 1 or die "Usage: $usage\n";

my $src_genome_id = shift;

my $id_server = Bio::KBase::IDServer::Client->new('http://bio-data-1.mcs.anl.gov/services/idserver');
my $cs = Bio::KBase::CDMI::Client->new('http://bio-data-1.mcs.anl.gov/services/cdmi_api');

my $out_fh;
if ($output_file)
{
    open($out_fh, ">", $output_file) or die "Cannot open $output_file: $!";
}
else
{
    $out_fh = \*STDOUT;
}

#my $genome_idx = $id_server->allocate_id_range('kb|g', 1);
#my $genome_id = "kb|g.$genome_idx";
#my $genome_id = "kb|g.22186";
my $genome_id = $src_genome_id;

my $ginfo = $cs->genomes_to_genome_data([$src_genome_id]);
$ginfo = $ginfo->{$src_genome_id};

my $domain = $ginfo->{taxonomy};
$domain =~ s/;.*$//;

my $contigs = [];
my $features = [];
my $genome = {
    id 		    => $genome_id,
    scientific_name => $ginfo->{scientific_name},
    domain 	    => $domain,
    genetic_code    => $ginfo->{genetic_code},
    contigs         => $contigs,
    features        => $features,
};

my $fres = $cs->get_relationship_IsOwnerOf([$src_genome_id], [], [], ['id', 'feature_type', 'source_id', 'function']);
my $fitems = [ map { $_->[2] } @$fres ];

my $fids = [ map { $_->{id} } @$fitems ];
my $trans = $cs->get_relationship_Produces($fids, [], [], ['id', 'sequence']);
my $f2trans = { map { $_->[1]->{from_link} => $_->[2] } @$trans };

my $locs = $cs->get_relationship_IsLocatedIn($fids, [], ['begin', 'dir', 'len', 'ordinal', 'to_link'], []);
my $f2loc = {};
for my $l (@$locs)
{
    my $id = $l->[1]->{from_link};
    my $item = $l->[1];
    
    $f2loc->{$id}->[$item->{ordinal}] = $item;
}

my $c = $cs->genomes_to_contigs([$src_genome_id]);
my $src_contigs = $c->{$src_genome_id};
my $cseqs = $cs->contigs_to_sequences($src_contigs);
my $cprefix = "$genome_id.c";
#my $cstart = $id_server->allocate_id_range($cprefix, scalar @$src_contigs) + 0;
my %cmap;
for my $ctg (@$src_contigs)
{
#    my $nctg = "$cprefix.$cstart";
#    $cstart++;
#    push(@$contigs, { id => $nctg, dna => $cseqs->{$ctg} });
#    $cmap{$ctg} = $nctg;
    $cmap{$ctg} = $ctg;
}

my %items;
for my $item (@$fitems)
{
    push(@{$items{$item->{feature_type}}}, $item);
}

for my $type (keys %items)
{
    my $fids = $items{$type};
    my $n = @$fids;
    my $prefix = "$genome_id.$type";
#    my $start = $id_server->allocate_id_range($prefix, $n) + 0;

    for my $item (sort { my $afid = $a->{id};
			 my $bfid = $b->{id};
			 my $aloc = $f2loc->{$afid}->[0];
			 my $bloc = $f2loc->{$bfid}->[0];
			 $aloc->{to_link} cmp $bloc->{to_link} or
			     $aloc->{begin} <=> $bloc->{begin} }
		  @$fids)
    {
	my $feature = {};
	my $fid = $item->{id};
#	my $nfid = "$prefix.$start";
	my $nfid = $fid;
#	$start++;
	my $locs = $f2loc->{$fid};

	my @nloc;
	for my $loc (@$locs)
	{
	    my $nctg = $cmap{$loc->{to_link}};

	    push(@nloc, [$nctg, $loc->{begin}, $loc->{dir}, $loc->{len}]);
	}

	if (exists($f2trans->{$fid}))
	{
	    my $trans = $f2trans->{$fid};
	    $feature->{protein_translation} = $trans->{sequence};
	}
	$feature->{location} = [@nloc];
	$feature->{function} = $item->{function};
	$feature->{id} = $nfid;
	$feature->{type} = $type;
	$feature->{annotations} = [];
	$feature->{aliases} = [];
	push(@$features, $feature);
    }
}
    
my $jdumper = JSON::XS->new;
$jdumper->pretty(1);
# Note - I added utf8 because without it the JSON reader on the other end chokes...
print $out_fh $jdumper->utf8->encode($genome);
close($out_fh);
