This repo contains scripts to compute (roughly) some confidence scores
for particular sets of roles that each protein in a query organism could play.

The wrapper script for computing probabilistic annotations from a workspace is probanno-driver:

$ probanno-driver [genome_id] [genome_workspace] [probanno_id] [probanno_workspace]

You must have at least r permission in genome_workspace and at least w permission in
probanno_workspace to successfully run this script.

This script takes a couple hours to run for an average-sized bacterial genome, so
it should eventually go into the queuing system Chris has set up for gapfill and other
long-running problems.

The package expects to be able to place its results and intermediate files
in a standard location: 

/kb/deployment/data/probabilistic_annotation

Results specific to a genome are saved in:

/kb/deployment/data/probabilistic_annotation/[genome_id]/

The intermediate files in /kb/deployment/data/probabilistic_annotation
take a nontrivial amount of time to generate and many of them are the same for 
any given query genome, so I recommend keeping the files around so that
next time the scripts are run it is much faster. They can be updated whenever KBase is updated
by re-running the ExtractorDriver.py script.

The probability algorithm is not yet finalized but it should give reasonable results
that can be used to test other functions that can depend on it. The results will also get
better as the database loads become more complete.
