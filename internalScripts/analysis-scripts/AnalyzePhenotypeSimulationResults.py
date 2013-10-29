#!/usr/bin/python

import json
import optparse
import sys

from biokbase.workspaceService.Client import *
from biokbase.fbaModelServices.Client import *

usage = "%prog -p [Phenotype ID] -w [Workspace] -a [Auth string] (options)"
description = '''

Given a JSON object calculate the false positive, negative rates. Optionally omit everything that isn't
in the specified media.

'''

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-a", "--auth", help="Auth token (required)", action="store", type="str", dest="auth", default=None)
parser.add_option("-w", "--phenows", help="Workspace for PhenotypeSimulationSet object (required)", action="store", type="str", dest="phenows", default=None)
parser.add_option("-p", "--pheno", help="ID for PhenotypeSimulationSet object (required)", action="store", type="str", dest="pheno", default=None)
parser.add_option("-m", "--media", help="(optional) Limit analysis to only this media condition", action="store", type="str", dest="media", default=None)
# Note - the model is already linked to the phenotypeSimulationSet but the rxnprobs is not
parser.add_option("-r", "--probrxn", help="(optional) Specify a RxnProbs object to get the probabilities for reactions associed with false-positive \ negative genes. [only applicable to knockout data]",
                  action="store", type="str", dest="probrxn", default=None)
parser.add_option("-x", "--probrxnws", help="(optional) RxnProbs workspace. Assumed to be the same as phenows if not specified.",
                  action="store", type="str", dest="probrxnws", default=None)
parser.add_option("-u", "--url", help="URL", action="store", type="str", dest="url", default="http://bio-data-1.mcs.anl.gov/services/fba")
(options, args) = parser.parse_args()

if options.auth is None or options.phenows is None or options.pheno is None:
    sys.stderr.write("ERROR: -a (auth), -w (phenows) and -p (pheno) are required arguments.\n")
    exit(2)

if options.probrxn is not None and options.probrxnws is None:
    options.probrxnws = options.phenows

# Get simulation data
wsClient = workspaceService("http://kbase.us/services/workspace/")
phenoobj = wsClient.get_object( { "workspace" : options.phenows,
                                  "type" : "PhenotypeSimulationSet",
                                  "id"   : options.pheno,
                                  "auth" : options.auth
                                  })

sims = phenoobj["data"]["phenotypeSimulations"]

modeldict = {}
wasgapfilled = {}
rxnprobdict = {}
if options.probrxn is not None:
    rxnprobobj = wsClient.get_object( { "workspace" : options.probrxnws,
                                        "type"      : "RxnProbs",
                                        "id"        : options.probrxn,
                                        "auth"      : options.auth
                                        })
    # If we specified a rxnProbs object we need the model object as well to use to grab the reaction IDs associated with each gene
    modelid = phenoobj["data"]["model"]
    modelws = phenoobj["data"]["model_workspace"]
    try:
        fbaClient = fbaModelServices(options.url)
        models = fbaClient.get_models( { "models" : [ modelid ],
                                         "workspaces" : [ modelws ],
                                         "auth"   : options.auth
                                         })
    except:
        sys.stderr.write("Unable to get model %s out of workspace %s - which is necessary to do rxnprobs-based analyses\n" %(modelid, modelws) )
        exit(0)
    model = models[0]
    # What we really need from these two objects are two dictionaries:
    # MODEL: gene -> reaction(s), gapfilledornot
    # RXNPROBS: reaction -> rxnprob_array ( [ rxnid, probability, diagnostic, complexes, GPR ] )
    for reaction in model["reactions"]:
        features = reaction["features"]
        for feature in features:
            if feature in modeldict:
                modeldict[feature].append( [ reaction["reaction"], reaction["gapfilled"] ] )
            else:
                modeldict[feature] = [ [ reaction["reaction"], reaction["gapfilled"] ] ]
                pass
            pass
        pass
    for rxnprob in rxnprobobj["data"]["reaction_probabilities"]:
        rxnprobdict[rxnprob[0]] = rxnprob
        pass    
    pass

# Set up counters...
CN = 0
CP = 0
FN = 0
FP = 0

if options.media is None:
    print "For phenotypedata %s and media ALL: " %(options.pheno)
else:
    print "For phenotypedata %s and media %s: " %(options.pheno, options.media)
    pass

# Count the CP\CN\FP\FN
# If a rxnprobs object is specified, print the rxnprobs data associated with each gene...
for sim in sims:
    media = sim[0][1]
    sim_type = sim[3]
    genes = sim[0][0]
    if options.media is not None:
        if options.media != media:
            continue
    if len(rxnprobdict.keys()) > 0:
        for gene in genes:
            # Some genes are not in the model
            if gene not in modeldict:
                continue
            for rxnpair in modeldict[gene]:
                # Some reactions in the model have no probabilities associated (we give them 0)
                rxn = rxnpair[0]
                isgapfilled = str(rxnpair[1])
                if rxn in rxnprobdict:
                    rxnprobobj = rxnprobdict[rxn]
                else:
                    rxnprobobj = [ rxn, 0, "NOCOMPLEX", "", "" ]
                # Gene, WasGapfilled, reaction, probability, diagnostics...
                print "\t".join( [ gene, isgapfilled ] + [ str(s) for s in rxnprobobj ] + [ sim_type ] )

    if sim_type == "CN":
        CN +=1
    elif sim_type == "CP":
        CP += 1
    elif sim_type == "FN":
        FN += 1
    elif sim_type == "FP":
        FP += 1
    pass

print "Total CN (correct negative): %d" %(CN)
print "Total CP (correct positive): %d" %(CP)
print "Total FN (false negative): %d" %(FN)
print "Total FP (false positive): %d" %(FP)
