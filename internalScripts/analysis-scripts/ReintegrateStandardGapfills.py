#!/usr/bin/python

from biokbase.fbaModelServices.Client import *
from biokbase.workspaceService.Client import *
from biokbase.probabilistic_annotation.Client import _read_inifile
import optparse
import os
import sys

usage = "%prog -m [gapfilled_model] -w [model_workspace] (other options)"
description = """
Take already-integrated gapfill solutions from a model that was NOT integrated with a rxnprobs object
and re-integrate it with a rxnprobs object (same reactions, but add gene associations to new reactions).
This should be done in order to enable a knockout lethality analysis of genes that would be added upon
post-processing a network-based gap fill solution. Solution 0 is assumed to be the one you want to integrate.
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--modelid", help="Model ID", action="store", type="str", dest="modelid", default=None)
parser.add_option("-w", "--ws", help="Workspace for model and for RxnProbs object if specified...", action="store", type="str", dest="ws", default=None)
parser.add_option("-u", "--url", help="URL for FBA model services", action="store", type="str", dest="url", default="http://bio-data-1.mcs.anl.gov/services/fba")
parser.add_option("-r", "--rxnprobsid", help="Rxnprobs object ID", action="store", type="str", dest="rxnprobsid", default=None)
(options, args) = parser.parse_args()

if options.modelid is None:
    raise IOError("modelid is required input")

if options.ws is None:
    raise IOError("Workspace is required input")

if options.rxnprobsid is None:
    raise IOError("Rxnprobsid is a required input")

authdata = _read_inifile()
token = authdata['token']
fbaClient = fbaModelServices(options.url)
wsClient = workspaceService("http://kbase.us/services/workspace/")

### Get the model object
#
models = fbaClient.get_models( { "models" : [ options.modelid ],
                                 "workspaces" : [ options.ws ],
                                 "auth"   : token
                                 })

gapfills = models[0]["integrated_gapfillings"]
if len(gapfills) < 1:
    raise IOError("ERROR: No integrated gapfillings found. Are you sure you are using the integrated model?")
print gapfills

solutionString = ";".join( [ s[1] + ".solution.0" for s in gapfills ] )

cmd = "kbfba-integratesolution -e -r %s -x %s -w %s -f \"%s\" -i %s %s" %(options.rxnprobsid, options.ws, options.ws, solutionString, options.modelid + ".REINTEGRATED", options.modelid)
print cmd
os.system(cmd)
