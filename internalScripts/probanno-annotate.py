#!/usr/bin/python

import optparse, sys
from biokbase.probabilistic_annotation.Client import *

usage="%prog [options]"
description="""Build a probabilistic annotation object"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--genome", help="ID of genome object", action="store", dest="genome", default=None)
parser.add_option("-e", "--genomews", help="ID of genome workspace", action="store", dest="genome_workspace", default=None)
parser.add_option("-p", "--probanno", help="ID of probanno object)", action="store", dest="probanno", default=None)
parser.add_option("-w", "--probannows", help="ID of probanno workspace", action="store", dest="probanno_workspace", default=None)
parser.add_option("-o", "--overwrite", help="Overwrite existing probanno object", action="store", dest="overwrite", default=False)
parser.add_option("-d", "--debug", help="Save generated data for debug purposes", action="store", dest="debug", default=False)
parser.add_option("-v", "--verbose", help="Print verbose output", action="store", dest="verbose", default=False)
parser.add_option("-a", "--auth", help="Auth token", action="store", dest="auth", default=None)
parser.add_option("-u", "--url", help="Server url", action="store", dest="url", default="http://localhost:7073")
(options, args) = parser.parse_args()

params = { "genome": options.genome, "genome_workspace": options.genome_workspace,
           "probanno": options.probanno, "probanno_workspace": options.probanno_workspace,
           "overwrite": options.overwrite, "debug": options.debug, "auth": options.auth }

client = ProbabilisticAnnotation(options.url)
jobid = client.annotate(params)
print("Job %s submitted" %(jobid))
exit(0)
