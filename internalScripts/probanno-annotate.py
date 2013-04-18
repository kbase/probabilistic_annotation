#!/usr/bin/python

import optparse, sys
from biokbase.probabilistic_annotation.MyImpl import annotate

usage="%prog [options]"
description="""Main driver to get data needed out of the KBase and store it locally.
Data will be stored in a local database autoReconInfo
All of this data is QUERY INDEPENDENT. It should all be the same
for any organism for which you want to do a reconstruction..."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--genome", help="ID of genome object", action="store", dest="genome", default=None)
parser.add_option("-e", "--genomews", help="ID of genome workspace", action="store", dest="genome_workspace", default=None)
parser.add_option("-p", "--probanno", help="ID of probanno object)", action="store", dest="probanno", default=None)
parser.add_option("-w", "--probannows", help="ID of probanno workspace", action="store", dest="probanno_workspace", default=None)
parser.add_option("-o", "--overwrite", help="Overwrite existing probanno object", action="store", dest="overwrite", default=False)
parser.add_option("-d", "--debug", help="Save generated data for debug purposes", action="store", dest="debug", default=False)
parser.add_option("-v", "--verbose", help="Print verbose output", action="store", dest="verbose", default=False)
parser.add_option("-a", "--auth", help="Auth token", action="store", dest="auth", default=None)
(options, args) = parser.parse_args()

success = annotate(options)
if success:
    exit(0)
exit(1)
