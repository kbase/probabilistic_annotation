#!/usr/bin/python

from biokbase.workspaceService.Client import *
from biokbase.probabilistic_annotation.Client import _read_inifile

import optparse
import sys

class SensitivityAnalyzer:
    def __init__(self, args):
        _authdata = _read_inifile()
        self.wsClient = workspaceService("http://kbase.us/services/workspace/")
        self.token = _authdata['token']
        self.args = args
    def _getObject(self, objType, objId, objWs):
        getObjectArgs = dict()
        getObjectArgs['type'] = objType
        getObjectArgs['id'] = objId
        getObjectArgs['workspace'] = objWs
        getObjectArgs['auth'] = self.token
        obj = self.wsClient.get_object(getObjectArgs)
        return obj

    def _getRxnSensitivityScores(self, RxnSensitivity):
        ''' _getRxnSensitivityScores - Calculate Reaction Sensitivity scores.

        This function should be run on a RxnSensitivity object that is generated with delete_reactions=FALSE.

        Reaction Sensitivity scores are computed as follows:
        Gapfill Reaction --> Inactivated reactions a and b

        if inactive reaction a was inactivated by 3 gapfill reactions and b by 2, the total cost is
        1/3 + 1/2 = 5/6

        '''

        # How many things does each inactive reaction fill?
        inactiveRxnsNumGapfilled = dict()
        gapfilledReactions = set()
        for reaction in RxnSensitivity['data']['reactions']:
            gapfilledReactions.add(reaction['reaction'])

        for reaction in RxnSensitivity['data']['reactions']:
            for inactive_rxn in reaction['new_inactive_rxns']:
                if inactive_rxn in gapfilledReactions:
#                    sys.stderr.write("WARNING: Inactive reaction %s for gapfill reaction %s was also a gapfilled reacton\n" %(inactive_rxn, reaction['reaction']))
                    inactiveRxnsNumGapfilled[inactive_rxn] = 0
                    continue
                if inactive_rxn in inactiveRxnsNumGapfilled:
                    inactiveRxnsNumGapfilled[inactive_rxn] += 1
                else:
                    inactiveRxnsNumGapfilled[inactive_rxn] = 1
        # The total cost is sum(1/N) where N is the number of things each gapfill reaction inactiveted
        gapfillRxnToCost = dict()
        for reaction in RxnSensitivity['data']['reactions']:
            cost = 0
            for inactive_rxn in reaction['new_inactive_rxns']:
                if inactiveRxnsNumGapfilled[inactive_rxn] != 0:
                    cost += float(1)/float(inactiveRxnsNumGapfilled[inactive_rxn])
            gapfillRxnToCost[reaction['reaction']] = cost

        return gapfillRxnToCost
                    

    def calculateSensitivityStatistics(self):
        # Lets first get basic data from the RxnSensitivity object
        RxnSensitivity = self._getObject('RxnSensitivity', self.args.RxnSensitivity_id, self.args.RxnSensitivity_ws)

        RxnProbsDict = dict()
        if self.args.RxnProbs_id is not None:
            RxnProbs = self._getObject("RxnProbs", self.args.RxnProbs_id, self.args.RxnProbs_ws)
            for reaction in RxnProbs['data']['reaction_probabilities']:
                rxnid = reaction[0]
                probability = reaction[1]
                RxnProbsDict[rxnid] = float(probability)

        gapfillRxnToCost = self._getRxnSensitivityScores(RxnSensitivity)

        # Now for the fun part.
        Header = ( "Gapfill_reaction_id", "Gapfill_reaction_likelihood", "Deleted", "Inactive_reactions", "Inactive_rxn_cost", "Average_inactive_rxn_likelihood" )
        FinalOutput = []
        for reaction in RxnSensitivity['data']['reactions']:
            rxnid = reaction['reaction']
            if rxnid in RxnProbsDict:
                rxn_probability = RxnProbsDict[rxnid]
            else:
                rxn_probability = 0
            inactive_rxn_list = reaction['new_inactive_rxns']
            inactive_rxn_string = ";".join(inactive_rxn_list)
            inactive_rxn_cost = gapfillRxnToCost[rxnid]
            n_inactive_rxns = len(inactive_rxn_list)
            rxn_deleted = reaction['deleted']
            # How likely were the inactivated reactions to begin with?
            if inactive_rxn_cost == 0:
                FinalOutput.append( ( rxnid, rxn_probability, rxn_deleted, inactive_rxn_string, inactive_rxn_cost, None) )
                continue
            else:
                inactive_rxn_sum = 0
                for inactive_rxn in inactive_rxn_list:
                    if inactive_rxn in RxnProbsDict:
                        inactive_rxn_sum += RxnProbsDict[inactive_rxn]
                average_inactive_rxn_likelihood = inactive_rxn_sum/float(n_inactive_rxns)
                FinalOutput.append( (rxnid, rxn_probability, rxn_deleted, inactive_rxn_string, inactive_rxn_cost, average_inactive_rxn_likelihood) )
        return FinalOutput, Header
        

if __name__ == "__main__":
    usage = "%prog [options] RxnSensitivity RxnSensitivityWorkspace"
    description = """Gather statistics on reactions in a rxnsensitivity object"""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("-r", "--rxnprobs", help="Rxnprobs object - probabilities for reactions will be reported along with reaction sensitivt statistics",
                      action="store", type="str", dest="RxnProbs_id", default=None)
    parser.add_option("-w", "--rxnprobsws", help="Workspace for RxnProbs object (default - same as RxnSensitivity)",
                      action="store", type="str", dest="RxnProbs_ws", default=None)
    (options, args) = parser.parse_args()

    if len(args) < 2:
        print usage
        exit(1)

    options.RxnSensitivity_id = args[0]
    options.RxnSensitivity_ws = args[1]

    if options.RxnProbs_ws is None and options.RxnProbs_id is not None:
        options.RxnProbs_ws = options.RxnSensitivity_ws

    SensitivityObject = SensitivityAnalyzer(options)
    SensitivityAnalysisResults, Header = SensitivityObject.calculateSensitivityStatistics()

    print "\t".join(Header)
    for res in SensitivityAnalysisResults:
        print "\t".join( [ str(s) for s in res ] )
