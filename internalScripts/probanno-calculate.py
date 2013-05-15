#!/usr/bin/python

import optparse, sys
from biokbase.probabilistic_annotation.Client import *

usage="%prog [options]"
description="""Calculate reaction probabilities"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-p", "--probanno", help="ID of probanno object)", action="store", dest="probanno", default=None)
parser.add_option("-w", "--probannows", help="ID of probanno workspace", action="store", dest="probanno_workspace", default=None)
parser.add_option("-d", "--debug", help="Save generated data for debug purposes", action="store", dest="debug", default=False)
parser.add_option("-v", "--verbose", help="Print verbose output", action="store", dest="verbose", default=False)
parser.add_option("-a", "--auth", help="Authentication token", action="store", dest="auth", default=None)
parser.add_option("-u", "--url", help="Probanno Server url", action="store", dest="url", default="http://localhost:7073")
(options, args) = parser.parse_args()

params = { "probanno": options.probanno, "probanno_workspace": options.probanno_workspace,
           "debug": options.debug, "verbose": options.verbose, "auth": options.auth }

client = ProbabilisticAnnotation(options.url)
rxnprobs = client.calculate(params)

for prob in rxnprobs:
    print "\t".join([ str(s) for s in prob ] )

exit(0)
