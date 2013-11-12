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

        # Now for the fun part.
        Header = ( "Gapfill_reaction_id", "Gapfill_reaction_likelihood", "Deleted", "Inactive_reactions", "N_inactive", "Average_inactive_rxn_likelihood" )
        FinalOutput = []
        for reaction in RxnSensitivity['data']['reactions']:
            rxnid = reaction['reaction']
            if rxnid in RxnProbsDict:
                rxn_probability = RxnProbsDict[rxnid]
            else:
                rxn_probability = 0
            inactive_rxn_list = reaction['new_inactive_rxns']
            inactive_rxn_string = ";".join(inactive_rxn_list)
            n_inactive_rxns = len(inactive_rxn_list)
            rxn_deleted = reaction['deleted']
            # How likely were the inactivated reactions to begin with?
            if n_inactive_rxns == 0:
                FinalOutput.append( ( rxnid, rxn_probability, rxn_deleted, inactive_rxn_string, n_inactive_rxns, None) )
                continue
            else:
                inactive_rxn_sum = 0
                for inactive_rxn in inactive_rxn_list:
                    if inactive_rxn in RxnProbsDict:
                        inactive_rxn_sum += RxnProbsDict[inactive_rxn]
                average_inactive_rxn_likelihood = inactive_rxn_sum/float(n_inactive_rxns)
                FinalOutput.append( (rxnid, rxn_probability, rxn_deleted, inactive_rxn_string, n_inactive_rxns, average_inactive_rxn_likelihood) )
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
    print options

    SensitivityObject = SensitivityAnalyzer(options)
    SensitivityAnalysisResults, Header = SensitivityObject.calculateSensitivityStatistics()

    print "\t".join(Header)
    for res in SensitivityAnalysisResults:
        print "\t".join( [ str(s) for s in res ] )
