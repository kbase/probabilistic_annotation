This repo contains scripts to compute (roughly) some confidence scores
for particular sets of roles that each protein in a query organism could play.

There are several ways this can be run:

- Using the probanno-driver script with an extant KBase Genome Identifier (e.g. kb|g.0)
- Using the Perl API (probabilistic_annotationImpl.pm) functions which take either a
  genomeTO object or a KBase genome ID.

If you are running the probanno-driver script MAKE SURE you use the following syntax:

$ probanno-driver "kb|g.0"

The package expects to be able to place its results and intermediate files
in a standard location (by default,
/kb/deployment/data/probabilistic_annotation ). Results specific to a genome are saved in 
/kb/deployment/data/probabilistic_annotation/[genome_id]/ .

These intermediate files take a nontrivial amount of time to generate and many of them are the same for 
any given query genome, so I recommend keeping the files around so that
next time the scripts are run it is much faster. They can be updated whenever KBase is updated.

The probability algorithm is not yet finalized but it should give reasonable results
that can be used to test other functions that can depend on it. The results will also get
better as the database loads become more complete.
