This repo contains scripts to compute (roughly) some confidence scores
for particular sets of roles that each protein in a query organism could play.

This is a python service, it requires installation of uwsgi and simplejson. These should
be included in the runtime.

Scripts:

pa-annotate  : Calculate probabilistic annotations from a genome
pa-calculate : Convert probabilistic annoations into reaction probabilities
pa-gendata   : Generate basis data needed for calculating probabilistic annotations 
	       and reaction probabilities
pa-normalize : OUTDATED
pa-url       : Get the currently-set URL for the probanno service
