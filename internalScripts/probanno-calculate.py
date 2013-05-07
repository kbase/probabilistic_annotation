#!/usr/bin/python

import optparse, sys
#from biokbase.probabilistic_annotation.Impl import ProbabilisticAnnotation
from biokbase.probabilistic_annotation.Client import *

usage="%prog [options]"
description="""Calculate reaction probabilities"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--model", help="ID of model object", action="store", dest="model", default=None)
parser.add_option("-q", "--modelws", help="ID of genome workspace", action="store", dest="model_workspace", default=None)
parser.add_option("-p", "--probanno", help="ID of probanno object)", action="store", dest="probanno", default=None)
parser.add_option("-w", "--probannows", help="ID of probanno workspace", action="store", dest="probanno_workspace", default=None)
parser.add_option("-o", "--overwrite", help="Overwrite existing probanno object", action="store", dest="overwrite", default=False)
parser.add_option("-d", "--debug", help="Save generated data for debug purposes", action="store", dest="debug", default=False)
parser.add_option("-v", "--verbose", help="Print verbose output", action="store", dest="verbose", default=False)
parser.add_option("-a", "--auth", help="Authentication token", action="store", dest="auth", default=None)
parser.add_option("-u", "--url", help="Server url", action="store", dest="url", default="http://localhost:7073")
(options, args) = parser.parse_args()

params = { "model": options.model, "model_workspace": options.model_workspace,
           "probanno": options.probanno, "probanno_workspace": options.probanno_workspace,
           "overwrite": options.overwrite, "debug": options.debug, "verbose": options.verbose,
           "auth": options.auth }

#impl_ProbabilisticAnnotation = ProbabilisticAnnotation(None)
#success = impl_ProbabilisticAnnotation.calculate(params)

client = ProbabilisticAnnotation(options.url)
object_meta = client.calculate(params)
exit(0)
