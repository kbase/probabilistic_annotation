#!/usr/bin/python

from biokbase.fbaModelServices.Client import *
from biokbase.workspaceService.Client import *
from biokbase.probabilistic_annotation.Client import _read_inifile
import optparse
import re
import sys

usage = "%prog -m [gapfilled_model] -w [model_workspace] (other options)"
description = """Given a GAPFILLED model, pull out the Gapfill solution 
(including objetive values of all alternative solutions available in the gapfill run"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--modelid", help="Model ID", action="store", type="str", dest="modelid", default=None)
parser.add_option("-w", "--ws", help="Workspace for model and for RxnProbs object if specified...", action="store", type="str", dest="ws", default=None)
parser.add_option("-u", "--url", help="URL for FBA model services", action="store", type="str", dest="url", default="http://bio-data-1.mcs.anl.gov/services/fba")
parser.add_option("-r", "--rxnprobsid", help="Rxnprobs object ID (optional)", action="store", type="str", dest="rxnprobsid", default=None)
(options, args) = parser.parse_args()

if options.modelid is None:
    raise IOError("modelid is required input")

if options.ws is None:
    raise IOError("Workspace is required input")

authdata = _read_inifile()
token = authdata['token']
fbaClient = fbaModelServices(options.url)
wsClient = workspaceService("http://kbase.us/services/workspace/")

### Get the model object
#
#try:
models = fbaClient.get_models( { "models" : [ options.modelid ],
                                 "workspaces" : [ options.ws ],
                                 "auth"   : token
                                 })
#except:
#    raise IOError("ERROR: Getting model %s from workspace %s failed (most likely this means the model does not exist in that workspace)" %(options.modelid, options.ws))

# Get IDs for all the reactions and genes in the model (pre-gapfilling)
rxns_in_model = set()
genes_in_model = set()
integrated_reactions = set()
for reaction in models[0]["reactions"]:
    if str(reaction["gapfilled"]) != "1":
        rxns_in_model.add(reaction["reaction"])
        for feature in reaction["features"]:
            genes_in_model.add(feature)
    else:
        integrated_reactions.add(reaction["reaction"])

### Get the gapfill UUID integrated into the model
# Assuming here that we only care about the first gapfill run
# If none are integrated, look for an unintegrated gapfill UUIC
gapfills = models[0]["integrated_gapfillings"]
if len(gapfills) < 1:
    sys.stderr.write("Integrated gapfillings not found. Trying unintegrated...\n")
    gapfills = models[0]["unintegrated_gapfillings"]
    if len(gapfills) < 1:
        raise IOError("Specified gapfill %s does not have any integrated or unintegrated gapfillings!" %(options.modelid))

# If a RxnProbs object is specified, use it to get the probabilities, GPRs and genes associated with
# each reaction
rxnToProbability = {}
rxnToGapfillGpr = {}
if options.rxnprobsid is not None:
    rxnprobs_object = wsClient.get_object( { "workspace"     : options.ws,
                                             "type"          : "RxnProbs",
                                             "id"            : options.rxnprobsid,
                                             "auth"          : token
                                             })
    for rxnprob in rxnprobs_object["data"]["reaction_probabilities"]:
        rxnToProbability[rxnprob[0]] = rxnprob[1]
        rxnToGapfillGpr[rxnprob[0]] = rxnprob[4]
        pass
    pass

# KBase
geneFinder = re.compile("kb\|g\.\d+\.peg\.\d+")
# SEED
geneFinder2 = re.compile("fig\|\d+\.\d+\.peg\.\d+")

for gapfill in gapfills:
    gapfill_uuid = gapfill[1]

    # Get the gapfill formulation object from the workspace (so that references to data are saved)
    # Note - FBA isn't technically the correct type but with getting objects
    # by reference the type doesn't matter...
    gapfill_object_workspace = wsClient.get_object( { "workspace"  : "NO_WORKSPACE",
                                                       "type"      : "FBA",
                                                       "id"        : gapfill_uuid,
                                                       "auth"      : token
                                                       })

    # We need to reach in and grab the FBA solution object. It contains the text of the ProblemReport.txt
    # that we need to parse to get the objective values out.
    try:
        fba_formulation_uuid = gapfill_object_workspace["data"]["fbaFormulation_uuid"]
    except KeyError:
        sys.stderr.write("WARNING: Gapfill with UUID %s had no FBA result attached to it...(in model %s)\n" %(gapfill_uuid, options.modelid))
        continue

    fba_formulation_object = wsClient.get_object( { "workspace" : "NO_WORKSPACE",
                                                    "type"      : "FBA",
                                                    "id"        : fba_formulation_uuid,
                                                    "auth"      : token
                                                    })

    # The problemreport.txt file contains the objective values...
    # But old versions of gapfill didn't return this file as part of the object.
    try:
        problem_report_text = fba_formulation_object["data"]["fbaResults"][0]["outputfiles"]["ProblemReport.txt"]
    except KeyError:
        sys.stderr.write("WARNING: No ProblemReport.txt was attached to FBA Results object %s in model %s\n" %(fba_formulation_uuid, options.modelid))
        continue

    #####################################################
    # PARSE ProblemReport.txt to get gapfill solutions  #
    #####################################################

    # For iterative gapfilling there is one problem report string for each target reaction that had a successful gapfill.
    # The first row, though, is a header row, so we just start at the second row.
    for problem_string in problem_report_text[1:]:
        # It has lots of stuff that is semicolon-delimited and the solution we want is the third column.
        spl_problem_string = problem_string.split(";")
        solution_string = spl_problem_string[2]
        # And it is formatted as follows... Recursive MILP [objective]:[+-][rxnid],[+-][rxnid],...| ( repeating).
        # So lets strip off the Recursive MILP...
        solution_string = solution_string.replace("Recursive MILP ", "")
        # Lets also remove any rogue spaces just in case.
        solution_string = solution_string.replace(" ", "")
        # Now we have pipe-delmited lists of solutions.
        solution_list = solution_string.split("|")
        solution_rxns_to_objective = {}
        for sln in solution_list:
            # Separate the objective value from the list of reactions.
            obj_list = sln.split(":")
            try:
                objective = float(obj_list[0])
            except ValueError:
                # The ProblemReport.txt can report empty solutions too (they just look like ")" or something). We don't want those.
                continue
            rxn_list_string = obj_list[1]
            # Remove + and - signs from the reactions.
            rxn_list_string = rxn_list_string.replace("+","").replace("-", "")
            # Sort the list so that we dont have to worry about the reactions
            # being in a different order from what we get below.
            rxn_list = sorted(rxn_list_string.split(","))
            solution_rxns_to_objective[tuple(rxn_list)] = objective

    ### Grab the gapfill object associated with that UUID using the FBA object
    # We are mostly using this as a sanity check to make sure that the reactions
    # match up between solutions gotten above and those returned by the FBA model service.
    #
    # It's likely we could just use what we have from the ProblemReport.txt instead.
    gapfill_objects_fba = fbaClient.get_gapfills( { "workspaces" : [ "NO_WORKSPACE" ],
                                                    "gapfills"   : [ gapfill_uuid ],
                                                    "auth"       : token
                                                    })

    ####################################################
    #    Iterate over gapfill solutions and generate   #
    #    a table of results                            #
    ####################################################
    gapfill_solutions = gapfill_objects_fba[0]["solutions"]
    ii = 0
    print "\t".join( [ "Sln_idx", "rxnid", "objective", "gapfill_uuid", "rxn_likelihood", "is_reversibility_change", "probanno_based_GPR", "New_genes", "Number_of_new_genes", "Rxn_integrated_and_not_deleted" ] )
    for solution in gapfill_solutions:
        reactionAdds = solution["reactionAdditions"]
        if len(reactionAdds) < 1:
            sys.stderr.write("WARNING: Solution had no reactions in it for model %s\n" %(options.modelid) )
            continue

        ### Get the reactions in that gapfill object.
        #
        rxnids = set()
        for addedrxn in reactionAdds:
            rxnids.add(addedrxn[0])

        rxnids = tuple(sorted(rxnids))
        if rxnids not in solution_rxns_to_objective:
#            raise IOError("Solution mismatch between ProblemReport.txt and the Gapfill object from kbfba-getgapfills in model %s" %(options.modelid))
            sys.stderr.write("WARNING: Solution mismatch between ProblemReport.txt and the Gapfill object from kbfba-getgapfills in model %s (this is only OK if you're analyzing complete gapfill)\n" %(options.modelid))

        for rxnid in rxnids:
            # Is the gapfilled reaction still in the model?
            if rxnid not in rxns_in_model and rxnid not in integrated_reactions:
                reactionInModel = "False"
            else:
                reactionInModel = "True"
                pass

            # Not all reactions are in the RxnProbs object.
            if rxnid in rxnToProbability:
                p = str(rxnToProbability[rxnid])
            else:
                p = "NO_PROBABILITY"
                pass

            # If the reaction was in the original model and is still in the gapfill solution, that means it was a reversibility change.
            if rxnid in rxns_in_model:
                revchange = "1"
            else:
                revchange = "0"
                pass

            # Get the GPR associated with the gapfilled reaction and figure out if the genes added are actually new genes or not.
            newgenes = []
            if rxnid in rxnToGapfillGpr:
                gpr = rxnToGapfillGpr[rxnid]
                for match in geneFinder.findall(gpr):
                    if match not in genes_in_model:
                        newgenes.append(match)
                        pass
                    pass
                for match in geneFinder2.findall(gpr):
                    if match not in genes_in_model:
                        newgenes.append(match)
                        pass
                    pass
            else:
                gpr = ""
                pass
            nnew = str(len(newgenes))
            newgenes = ";".join(newgenes)

            # This doesn't work for complete gapfilling. I just want to get a hack aroud it.
            try:
                print "\t".join( [ str(ii), rxnid, str(solution_rxns_to_objective[rxnids]), gapfill_uuid, p, revchange, gpr, newgenes, nnew, reactionInModel ] )
            except KeyError:
                print "\t".join( [ str(ii), rxnid, str(None), gapfill_uuid, p, revchange, gpr, newgenes, nnew, reactionInModel ] )
        
        ii += 1

