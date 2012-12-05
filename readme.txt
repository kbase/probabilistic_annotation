This repo contains scripts to compute (roughly) some confidence scores
for particular sets of roles that each protein in a query organism could play.

There are several ways this can be run:

- Using the test_impl.pl script with an extant KBase Genome Identifier (e.g. kb|g.0)
- Using the Perl API (probabilistic_annotationImpl.pm) functions which take either a
  genomeTO object or a KBase genome ID.

If you are running the test_impl.pl script MAKE SURE you use the following syntax:

$ perl test_impl.pl "kb|g.0"

and NOT:

$ ./test_impl.pl "kb|g.0"

(the latter does not work on the SEED machines becasue the install of perl is not in the standard location)

The scripts require generation of some intermediate files which are the same for
any organism that you want to run against. These intermediate files take a nontrivial
amount of time to generate, so I recommend running and then saving the files so that
next time the scripts are run it is much faster.

The scripts at the moment do some backhanded things but it will be cleaned somewhat up once
we get a good API routine for cs_to_genomes up and running. 

The probability algorithm is not yet finalized but it should give reasonable results
that can be used to test other functions that can depend on it. The results will also get
better as the database loads become more complete.
