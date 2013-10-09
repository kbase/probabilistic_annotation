#!/usr/bin/python

import optparse
from biokbase.fbaModelServices.Client import *
from biokbase.workspaceService.Client import *

usage = "%prog -w workspace -a auth -m [Model ID] (options)"
description = ''' Produce a (model, gapfill_solutionid) tab-delimited file for use when integrating gapfill solutions. '''
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-s", "--solution", help="Solution number to integrate, starting from 0 (D: 0. Use -1 to get ALL solutions, which is recommended for iterative gapfilling)", action="store", type="str", dest="solution", default="0")
parser.add_option("-a", "--auth", help="auth token", action="store", type="str", dest="auth", default=None)
parser.add_option("-w", "--workspace", help="workspace in which model is found", action="store", type="str", dest="ws", default=None)
parser.add_option("-m", "--modelid", help="Model ID to search for gapfills", action="store", type="str", dest="modelid", default=None)
parser.add_option("-u", "--url", help="URL", action="store", type="str", dest="url", default="http://bio-data-1.mcs.anl.gov/services/fba")
(options, args) = parser.parse_args()

if options.ws is None or options.auth is None or options.modelid is None:
    raise IOError("workspace, auth and model are required arguments.")

fbaClient = fbaModelServices(options.url)
wsClient = workspaceService("http://kbase.us/services/workspace/")

try:
    models = fbaClient.get_models( { "models" : [ options.modelid ],
                                 "workspaces" : [ options.ws ],
                                 "auth"   : options.auth
                                 })
except:
    raise IOError("ERROR: Getting model %s from workspace %s failed (most likely this means the model does not exist in that workspace)" %(options.modelid, options.ws))

### Get the gapfill UUID integrated into the model
# Assuming here that we only care about the first gapfill run
# If none are integrated, look for an unintegrated gapfill UUIC
gapfills = models[0]["unintegrated_gapfillings"]
if len(gapfills) < 1:
    raise IOError("Specified gapfill %s does not have any unintegrated gapfillings!" %(options.modelid))

for gapfill in gapfills:
    gapfill_uuid = gapfill[1]
    if options.solution == -1:
        gapfill = wsClient.get_object( { "auth" : options.auth,
                                         "workspace" : options.ws,
                                         "workspace" : "NO_WORKSPACE",
                                         "type" : "GapFill",
                                         "id"   : gapfill_uuid
                                         })
        for ii in range(len(gapfill["data"]["gapfillingSolutions"])):
            gapfill_solutionid = "%s.solution.%s" %(gapfill_uuid, ii)
            print "%s\t%s" %(options.modelid, gapfill_solutionid)
        pass
    else:
        gapfill_solutionid = "%s.solution.%s" %(gapfill_uuid, options.solution)
        print "%s\t%s" %(options.modelid, gapfill_solutionid)
