#!/usr/bin/python

from CDMI import CDMI_EntityAPI
from DataExtractor import getFieldFromEntity, getFieldFromRelationship
# For the CDMI URL and other constants
from PYTHON_GLOBALS import *

import os
import sys

try:
    import json
except ImportError:
    # Needed for old version of python such as those on the SEED machines.
    # Add the simplejson folder to your PYTHONPATH variable or this will fail.
    import simplejson as json


def makeReactionProbabilityObject(reaction_probability_file, biochemistry_json_file, output_file):
    '''Searches the biochemistry JSON for reaction UUIDs.
    A dictionary is created (using the alias table for 'modelSEED' from reaction UUID
    to modelSEED ID, and another from modelSEED IDs to KBase reaction IDs.

    If we managed to get a probability for that reaction, we print that (even if it is 0 - which
    means that the complex was defined but not found in the organism) along with the proposed GPR
    which is just a string.

    A new object ('Reaction probability object') is created that contains the Biochemistry UUID and a list
    of (reaction UUID, probability, GPR, string, flags) sets that can be added to a model object'''

    if os.path.exists(output_file):
        sys.stderr.write("Modified biochemistry JSON file %s already exists!\n" %(output_file))
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
    # to each of these.
    # 
    # In the meantime, we make a new object to store this information so that we don't have to modify
    # the biochemistry object and keep it static.
    rxnProbObject = {}
    rxnProbObject["biochemistry_uuid"] = resp["uuid"]
    rxnProbObject["reaction_probabilities"] = []

    rxnList = resp["reactions"]
    for ii in range(len(rxnList)):
        myuuid = rxnList[ii]["uuid"]
        myProbability = 0
        myComplexInfo = ""
        myGPR = ""
        # These flags indicate database issues that could bias the probabilities.
        #
        # If all of them are FALSE or if only reactionHasComplex is true, the probability will be 0
        # but it is due to missing data.
        #
        # If only allComplexesHaveRepresentativeRoles is false, that only means there is some missing data
        # but still enough to test presence of some subunits of some complexes.
        reactionHasComplex = False
        oneComplexHasRepresentativeRoles = False
        allComplexesHaveRepresentativeRoles = False
        probablyNotThere = False
        # If the database versions are consistent this should always be the case
        if myuuid in uuidToModelId:
            modelId = uuidToModelId[myuuid]
            if modelId in modelToKbase:
                kbaseId = modelToKbase[modelId]
                # This one is only the case if there are complexes associated with the reaction
                if kbaseId in kbaseIdToInfo:
                    reactionHasComplex = True
                    ''' There are three possibilities for each complex.
                    1: The roles attached to the complex had no representatives in our BLAST db (BAD)
                    2: The roles attached to the complex had representatives, but they were not found in
                    our BLAST search (OK)
                    3: The roles attached were all found with some probability in the BLAST search (OK)

                    Since there are multiple possibilities for complexes we need to decide when we should
                    call it OK and when we can't.

                    For now I will set separate flags for the occasion when one complex has roles and one doesn't
                    (and the calculated probability is for the complex with a probability)
                    and the occasion when NONE of the complexes have representeatives of their roles
                    (and the calculated probability is artificially 0)

                    PARTIAL cases are treated as "has representatives" for this purpose.
                    Therefore, if only allComplexesHaveRepresentativeRoles is false, that means
                    there is incomplete information, but at least one subunit of one complex attached
                    to the reaction had a representative that we could use to calculate a probability.

                    '''
                    # CPLX_FULL   [ok]
                    # CPLX_NOTTHERE [ok]
                    # CPLX_PARTIAL [sort of ok - treated as OK for this purpose]
                    # CPLX_NOREPS [bad]
                    myProbability = kbaseIdToInfo[kbaseId][0]
                    myGPR = kbaseIdToInfo[kbaseId][2]
                    myComplexInfo = kbaseIdToInfo[kbaseId][1]
                    if "CPLX_NOREPS" in myComplexInfo:
                        if "CPLX_FULL" in myComplexInfo or "CPLX_NOTTHERE" in myComplexInfo or "CPLX_PARTIAL" in myComplexInfo:
                            oneComplexHasRepresentativeRoles = True
                        else:
                            # No complexes have representative roles.
                            pass
                    else:
                        # All of them are either CPLX_FULL or CPLX_NOTTHERE
                        oneComplexHasRepresentativeRoles = True
                        allComplexesHaveRepresentativeRoles = True
                    if "CPLX_NOTTHERE" in myComplexInfo and not ("CPLX_FULL" in myComplexInfo) and not ("CPLX_PARTIAL" in myComplexInfo) and not ("CPLX_NOREPS" in myComplexInfo):
                        probablyNotThere = True

        oneReactionObject = {}
        oneReactionObject["reaction_uuid"] = myuuid
        oneReactionObject["probability"] = myProbability
        oneReactionObject["complexinfo"] = myComplexInfo
        oneReactionObject["GPR"]         = myGPR
        oneReactionObject["reactionHasComplex"]                   = reactionHasComplex
        oneReactionObject["oneComplexHasRepresentativeRoles"]     = oneComplexHasRepresentativeRoles
        oneReactionObject["allComplexesHaveRepresentativeRoles"]  = allComplexesHaveRepresentativeRoles
        oneReactionObject["probablyNotThere"] = probablyNotThere
        # This isn't strictly necessary but makes it easier for me to interpret the data
        oneReactionObject["name"] = rxnList[ii]["name"]
        rxnProbObject["reaction_probabilities"].append(oneReactionObject)

    json.dump(rxnProbObject, open(output_file, "w"), indent=4)

if __name__ == "__main__":
    import optparse
    usage = "%prog [rxnprobfile] [biochemistry_json_file] [output_file]"
    description = """Uses the biochemistry object to identify reaction UUIDs corresponding to KBase IDs for which probabilities
were computed, and then makes a new (model reaction probability) object containing the probability information."""
    parser = optparse.OptionParser(usage=usage, description=description)
    (options, args) = parser.parse_args()
    if len(args) < 3:
        sys.stderr.write("ERROR: Incorrect usage of function\n")
        sys.stderr.write("USAGE: %s\n" %(usage))
        exit(2)       
    makeReactionProbabilityObject(args[0], args[1], args[2])

#addRxnProbabilitiesToBiochemistryJson("kb|g.0/kb|g.0.rxnprobs", "default_biochemistry.json", "RESULT.json")
