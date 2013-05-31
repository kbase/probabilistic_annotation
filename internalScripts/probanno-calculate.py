#!/usr/bin/python

import optparse, sys
from biokbase.probabilistic_annotation.Client import *

usage="%prog [options]"
description="""Calculate reaction probabilities"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-p", "--probanno", help="ID of probanno object", action="store", dest="probanno", default=None)
parser.add_option("-w", "--probannows", help="ID of probanno workspace", action="store", dest="probanno_workspace", default=None)
parser.add_option("-d", "--debug", help="Save generated data for debug purposes", action="store", dest="debug", default=False)
parser.add_option("-v", "--verbose", help="Print verbose output", action="store", dest="verbose", default=False)
parser.add_option("-a", "--auth", help="Authentication token (required)", action="store", dest="auth", default=None)
parser.add_option("-t", "--templatemodel", help="Template model (D: None - get dicts from CDMI)", action="store", dest="template_model", default=None)
parser.add_option("-m", "--templatemodelws", help="Template model workspace(D: None - get dicts from CDMI)", action="store", dest="template_model_workspace", default=None)
parser.add_option("-u", "--url", help="Probanno Server url", action="store", dest="url", default="http://localhost:7073")
parser.add_option("-o", "--outputid", help="ID for output reaction probability object", action="store", dest="outputid", default=None)
parser.add_option("-x", "--outputws", help="Output reaction probability workspace", action="store", dest="outputws", default=None)
(options, args) = parser.parse_args()

params = { "probanno"                 : options.probanno, 
           "probanno_workspace"       : options.probanno_workspace,
           "debug"                    : options.debug, 
           "verbose"                  : options.verbose, 
           "auth"                     : options.auth,
           "template_model"           : options.template_model, 
           "template_model_workspace" : options.template_model_workspace,
           "rxnprobs"                 : options.outputid,
           "outputws"                 : options.outputws }

client = ProbabilisticAnnotation(options.url)
metadata = client.calculate(params)

print metadata

exit(0)
