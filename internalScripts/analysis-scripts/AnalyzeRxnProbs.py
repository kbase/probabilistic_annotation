#! /usr/bin/python

import sys
import os
import argparse
from biokbase.workspaceService.client import *


# Main script function
if __name__ == "__main__":

    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='AnalyzeRxnProbs')
    parser.add_argument('rxnprobs', help='id of RxnProbs object', action='store')
    parser.add_argument('workspace', help='id of workspace', action='store')
    parser.add_argument('-?', '--usage', help='show usage information', action='store_true', dest='usage')
    parser.add_argument('-b', '--bins', type=int, help='number of bins for probability values', action='store', dest='numBins', default=10)
    args = parser.parse_args()
    
    if args.usage:
        print parser.usage()
        exit(0)
        
    # Get an authorization token.
    # MBM Can remove this when type compiler is fixed -- JIRA issue 282
    try:
        fid = open(os.path.join(os.environ["HOME"], ".kbase_auth"), "r")
        token = fid.read()
        token.strip()
        fid.close()
    except:
        sys.stderr.write("Failed to get an authorization token.\n")
        exit(1)
    
    # Create a workspace client.
    wsClient = workspaceService('http://www.kbase.us/services/workspace')
    
    # Load the RxnProbs object.
    rxnprobs = wsClient.get_object( { "type": "RxnProbs", "id": args.rxnprobs, "workspace": args.workspace, "auth": token } )
    
    binList = []
    for i in range(0,args.numBins+1):
        binList.append(0)
        
    # Run the list of reactions and put into bins.
    for reaction in rxnprobs['data']['reaction_probabilities']:
        if reaction[1] == 0:
            binList[args.numBins] += 1
        else:
            bin = int(reaction[1] * args.numBins)
            binList[bin] += 1
        
    # Print bins in csv format.
    print '0,%d' %(binList[args.numBins])
    for i in range(0,args.numBins):
        print '%4.2f,%d' %(1.0/args.numBins*(i+1), binList[i])
    
    exit(0)
    