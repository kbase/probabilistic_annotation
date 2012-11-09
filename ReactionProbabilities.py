#!/usr/bin/python

from CDMI import CDMI_EntityAPI
from DataExtractor import getFieldFromEntity, getFieldFromRelationship
from PYTHON_GLOBALS import *

import json
import os
import sys

def addRxnProbabilitiesToBiochemistryJson(reaction_probability_file, biochemistry_json_file, output_file):
    '''Searches the biochemistry JSON for reaction UUIDs.
    A dictionary is created (using the alias table for 'modelSEED' from reaction UUID
    to modelSEED ID, and another from modelSEED IDs to KBase reaction IDs.

    If we managed to get a probability for that reaction, we print that (even if it is 0 - which
    means that the complex was defined but not found in the organism.) along with the proposed GPR
    which is just a string.

    The probability of the reaction is in the 'probability' field while the GPR is in the 'GPR' field in
    the modified biochemistry json file.

    If we did NOT calculate a probability for a particular reaction that means no complexes are defined
    for it and we print -5 to indicate that those have '0 probability' but due to database limitations
    rather than due to lack of genetic evidence...'''

    if os.path.exists(output_file):
        sys.stderr.write("Modified biochemistry JSON file %s already exists!\n")
        exit(2)

    # KBase ID --> (probability, complex_info, GPR)
    kbaseIdToInfo = {}
    for line in open(reaction_probability_file, "r"):
        spl = line.strip("\r\n").split("\t")
        kbaseIdToInfo[spl[0]] = ( spl[1], spl[3], spl[4] )

    # Model ID --> Kbase ID
    cdmi_entity = CDMI_EntityAPI(URL)
    rxniddict = cdmi_entity.all_entities_Reaction(MINN, COUNT, ["source_id"])
    kbaseIds = getFieldFromEntity(rxniddict, "id")
    modelIds = getFieldFromEntity(rxniddict, "source_id")
    modelToKbase = {}
    for ii in range(len(modelIds)):
        modelToKbase[modelIds[ii]] = kbaseIds[ii]

    # So begins the biochemistry-dependent parts.
    # Different biochemistries will (I think?) have different UUIDs
    # for all the reactions in them... but they will have links set up
    # to the model IDs. At least, I HOPE so.
    resp = json.load(open(biochemistry_json_file, "r"))
    
    # UUID --> Model ID
    aliasSetList = resp["aliasSets"]
    uuidToModelId = {}
    for aliasSet in aliasSetList:
        if aliasSet["source"] == "ModelSEED" and aliasSet["attribute"] == "reactions":
            aliasDict = aliasSet["aliases"]
            for k in aliasDict:
                # aliasDict is really a dict from reaction id to a LIST of UUIDs, implying that
                # it is possible that more than one UUID goes with the same reaction ID.
                # If that happens (WHY?????) then I'll just assign all of them the probability
                # of that reaction.
                for uuid in aliasDict[k]:
                    uuidToModelId[uuid] = k
            # We found the one we need, no need to go through the rest of them...
            break
    
    # Now we need to iterate over all of the reactions and add the appropriate probabilities
    # to each of these...
    print uuidToModelId

    rxnList = resp["reactions"]
    for ii in range(len(rxnList)):
        myuuid = rxnList[ii]["uuid"]
        # -5 means database issue \ not found.
        myProbability = -5
        myComplexInfo = ""
        myGPR = ""
        if myuuid in uuidToModelId:
            print "uuidOK: %s" %(myuuid)
            modelId = uuidToModelId[myuuid]
            if modelId in modelToKbase:
                print "modelIdOK: %s" %(modelId)
                kbaseId = modelToKbase[modelId]
                print kbaseId
                if kbaseId in kbaseIdToInfo:
                    print "kbaseIdOK: %s" %(kbaseId)
                    myProbability = kbaseIdToInfo[kbaseId][0]
                    myComplexInfo = kbaseIdToInfo[kbaseId][1]
                    myGPR = kbaseIdToInfo[kbaseId][2]
        resp["reactions"][ii]["probability"] = myProbability
        resp["reactions"][ii]["complexinfo"] = myComplexInfo
        resp["reactions"][ii]["GPR"]         = myGPR
    json.dump(resp, open(output_file, "w"), indent=4)

addRxnProbabilitiesToBiochemistryJson("kb|g.0/kb|g.0.rxnprobs", "default_biochemistry.json", "RESULT.json")
