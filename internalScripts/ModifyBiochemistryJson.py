#!/usr/bin/python

from biokbase.cdmi.client import CDMI_EntityAPI
from biokbase.probabilistic_annotation.DataExtractor import getFieldFromEntity, getFieldFromRelationship
# For the CDMI URL and other constants
from biokbase.probabilistic_annotation.PYTHON_GLOBALS import *

import os
import sys

try:
    import json
except ImportError:
    # Needed for old version of python such as those on the SEED machines.
    # Add the simplejson folder to your PYTHONPATH variable or this will fail.
    import simplejson as json


def addRxnProbabilitiesToBiochemistryJson(reaction_probability_file, biochemistry_json_file, output_file):
    '''Searches the biochemistry JSON for reaction UUIDs.
    A dictionary is created (using the alias table for 'modelSEED' from reaction UUID
    to modelSEED ID, and another from modelSEED IDs to KBase reaction IDs.

    If we managed to get a probability for that reaction, we print that (even if it is 0 - which
    means that the complex was defined but not found in the organism) along with the proposed GPR
    which is just a string.

    The probability of the reaction is in the 'probability' field while the GPR is in the 'GPR' field in
    the modified biochemistry json file.

    If we did NOT calculate a probability for a particular reaction that means no complexes are defined
    for it and we print -5 to indicate that those have '0 probability' but due to database limitations
    rather than due to lack of genetic evidence...'''

    if os.path.exists(output_file):
        sys.stderr.write("Modified biochemistry JSON file %s already exists!\n" %(output_file))
        exit(2)

    # KBase ID --> (probability, complex_info, GPR)
    kbaseIdToInfo = {}
    for line in open(reaction_probability_file, "r"):
        spl = line.strip("\r\n").split("\t")
        kbaseIdToInfo[spl[0]] = ( spl[1], spl[3], spl[4] )

    # Model ID --> Kbase ID
    cdmi_entity = CDMI_EntityAPI(CDMI_URL)
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

        resp["reactions"][ii]["probability"] = myProbability
        resp["reactions"][ii]["complexinfo"] = myComplexInfo
        resp["reactions"][ii]["GPR"]         = myGPR
        resp["reactions"][ii]["reactionHasComplex"]                   = reactionHasComplex
        resp["reactions"][ii]["oneComplexHasRepresentativeRoles"]     = oneComplexHasRepresentativeRoles
        resp["reactions"][ii]["allComplexesHaveRepresentativeRoles"]  = allComplexesHaveRepresentativeRoles
        
    json.dump(resp, open(output_file, "w"), indent=4)

if __name__ == "__main__":
    import optparse
    usage = "%prog [rxnprobfile] [biochemistry_json_file] [output_file]"
    description = "Add reaction probabilities, complex information, and putative GPR to all reactions in a biochemistry object"
    parser = optparse.OptionParser(usage=usage, description=description)
    (options, args) = parser.parse_args()
    if len(args) < 3:
        sys.stderr.write("ERROR: Incorrect usage of function\n")
        sys.stderr.write("USAGE: %s\n" %(usage))
        exit(2)       
    addRxnProbabilitiesToBiochemistryJson(args[0], args[1], args[2])

#addRxnProbabilitiesToBiochemistryJson("kb|g.0/kb|g.0.rxnprobs", "default_biochemistry.json", "RESULT.json")
