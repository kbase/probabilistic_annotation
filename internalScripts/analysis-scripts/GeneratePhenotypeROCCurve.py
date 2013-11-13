#!/usr/bin/python -u

import argparse
import sys
import re
import os
import traceback
from operator import itemgetter, attrgetter
from biokbase.probabilistic_annotation.Client import _read_inifile
from biokbase.workspaceService.Client import *

desc = '''
Given a JSON object calculate the false positive, negative rates. Optionally omit everything that isn't
in the specified media.
'''
separator ='///'
rolesToRxnsFileName = 'roles_to_reactions'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='GeneratePhenotypeROCCurve.py', description=desc)
    parser.add_argument('genome', help='ID of genome to analyze', action='store', default=None)
    parser.add_argument('model', help='ID of integrated gap filled model to analyze', action='store', default=None)
    parser.add_argument('phenosimset', help='ID of PhenotypeSimulationSet object to analyze', action='store', default=None)
    parser.add_argument('workspace', help='ID of workspace containing objects to analyze', action='store', default=None)
    parser.add_argument('--solution', help='index number of solution in integrated gap filled model', action='store', default='0')
    parser.add_argument('--phenosimsetws', help='workspace containing PhenotypeSimulationSet object (same as workspace if not specified)', action='store', default=None)
    parser.add_argument('--media', help='limit analysis to only this media condition', action='store', dest='media', default=None)
    parser.add_argument('--probanno', help='ID of ProbAnno object for genome', action='store', dest='probanno', default=None)
    parser.add_argument('--probannows', help='workspace containing ProbAnno object (same as workspace if not specified)', action='store', dest='probannows', default=None)
    parser.add_argument('--rxnprobs', help='ID of RxnProbs object for genome', action='store', dest='rxnprobs', default=None)
    parser.add_argument('--rxnprobsws', help='workspace containing RxnProbs object (same as workspace if not specified)', action='store', dest='rxnprobsws', default=None)
    parser.add_argument('--script-dir', help='path to directory with analysis scripts', action='store', dest='scriptDir', default='.')
    parser.add_argument('--fba-url', help='url for fba modeling service', action='store', dest='fbaurl', default='http://bio-data-1.mcs.anl.gov/services/fba')
    parser.add_argument('--ws-url', help='url for workspace service', action='store', dest='wsurl', default='http://www.kbase.us/services/workspace/')
    args = parser.parse_args()
    
    if args.probanno is None:
        args.probanno = args.genome + '.probanno'
    if args.probannows is None:
        args.probannows = args.workspace
    if args.rxnprobs is None:
        args.rxnprobs = args.genome + '.rxnprobs'
    if args.rxnprobsws is None:
        args.rxnprobsws = args.workspace
    if args.phenosimsetws is None:
        args.phenosimsetws = args.workspace
        
    # Create the clients for the services we need.
    authdata = _read_inifile()
    token = authdata['token']
    wsClient = workspaceService(args.wsurl)
    
    # Build roles to reactions dictionary.
    step = 1
    print "+++ Step %d: Build reactions to roles dictionary +++" %(step)
    if not os.path.exists(rolesToRxnsFileName):
        os.system("all_roles_used_in_models | roles_to_complexes | get_relationship_HasStep -to id >%s" %(rolesToRxnsFileName))
    rolesToReactions = dict()
    reactionsToRoles = dict()
    # Each line of the file has four fields: (1) role, (2) hypothetical flag, (3) complex id, (4) reaction id
    rolesToRxnsFile = open(rolesToRxnsFileName, 'r')
    for line in rolesToRxnsFile:
        fields = line.strip('\r\n').split('\t')
        reaction = fields[3]
        if reaction not in reactionsToRoles:
            reactionsToRoles[reaction] = list()
        reactionsToRoles[reaction].append(fields[0])
    rolesToRxnsFile.close()
    print "  %d reactions in reactions to roles dictionary" %(len(reactionsToRoles))

    # Analyze the gap fill results for the specified model.
    step += 1
    print "+++ Step %d: Run AnalyzeGapfillResults.py for model '%s/%s' +++" %(step, args.workspace, args.model)
    rxnprobs = args.rxnprobs
    resultsFileName = args.model + '.results'
    os.system("%s/AnalyzeGapfillResults.py -m %s -w %s --rxnprobs %s --url %s >%s" \
              %(args.scriptDir, args.model, args.workspace, rxnprobs, args.fbaurl, resultsFileName))
    genesToReactions = dict()
    resultsFile = open(resultsFileName, 'r')
    resultsFile.readline() # Throw away header line
    for line in resultsFile:
        fields = line.strip('\r\n').split('\t')
        if fields[0] == args.solution:
            if fields[5] == '0': # Only keep reactions that are not a reversibility change
                if fields[6]: # Only keep reactions that have a GPR
                    rxnid = re.sub(r'rxn0*(\d+)', r'kb|rxn.\1', fields[1])
                    geneList = re.findall('fig\|\d+\.\d+\.peg\.\d+', fields[6])
                    for index in range(len(geneList)):
                        if geneList[index] not in genesToReactions:
                            genesToReactions[geneList[index]] = dict()
                        genesToReactions[geneList[index]][rxnid] = 0.0
    resultsFile.close()
    print "  %d genes in genes to reactions dictionary" %(len(genesToReactions))

    # Get the ProbAnno object from the workspace.
    step += 1
    probanno = args.probanno
    print "+++ Step %d: Get ProbAnno object '%s/%s'" %(step, args.probannows, probanno)
    paObject = wsClient.get_object( { 'id': probanno, 'workspace': args.probannows, 'type': 'ProbAnno', 'auth': token } )
    probAnno = paObject['data']
    print "  %d genes in ProbAnno roleset probabilities dictionary" %(len(probAnno['roleset_probabilities']))
    
    # Need to go through rolesets dictionary
    # If an entry has more than one role, split into parts
    # Then run the array and combine any duplicate roles by adding probs
    # Parse the roleset probabilities from the ProbAnno object.
    step += 1
    print "+++ Step %d: Parse rolesets into roles and adjust probabilities for duplicates +++" %(step)
    rolesetProbabilities = dict()
    for gene in probAnno['roleset_probabilities']:
        geneRoleList = probAnno['roleset_probabilities'][gene]
        geneRoleDict = dict()
        for index in range(len(geneRoleList)):
            prob = geneRoleList[index][1] # Probability for this roleset
            # Split multiple roles in roleset for this gene
            roleList = geneRoleList[index][0].split(separator)
            # If role occurs more than once, add up the probabilities
            for j in range(len(roleList)):
                if roleList[j] in geneRoleDict:
                    geneRoleDict[roleList[j]] += prob
                else:
                    geneRoleDict[roleList[j]] = prob
        rolesetProbabilities[gene] = geneRoleDict
    print "  %d genes in parsed roleset probabilities dictionary" %(len(rolesetProbabilities))                
    
    # for each reaction in the reactions dictionary, find the roles in the rolesToReactions dictionary
    # then find the roles in the probanno object for the gene
    step += 1
    print "+++ Step %d: Find maximum probability of reaction given gene +++" %(step)
    probsFile = open(args.genome+'.probs', 'w')
    numProbs = 0
    for gene in genesToReactions:
        if gene in rolesetProbabilities:
            geneRoleDict = rolesetProbabilities[gene]
            for reaction in genesToReactions[gene]:
                if reaction not in reactionsToRoles:
                    print 'Why is reaction %s not in the reactionToRoles dictionary?' %(reaction)
                roleList = reactionsToRoles[reaction]
                for index in range(0,len(roleList)):
                    for role in geneRoleDict:
                        if role == roleList[index]:
                            probsFile.write('P(%s | %s) = %f\n' %(reaction, gene, geneRoleDict[role]))
                            genesToReactions[gene][reaction] = max(geneRoleDict[role], genesToReactions[gene][reaction])
                            numProbs += 1
        else:
            print 'Gene %s not found in ProbAnno object' %(gene)
    probsFile.close()
    print "  %d reaction probabilities set in genesToReactions dictionary" %(numProbs)
        
    # Get the PhenotypeSimulationSet object.
    step += 1
    print "+++ Step %d: Get phenotype simulation set %s/%s +++" %(step, args.phenosimsetws, args.phenosimset)
    pssObject = wsClient.get_object( { 'id': args.phenosimset, 'workspace': args.phenosimsetws, 'type': 'PhenotypeSimulationSet', 'auth': token } )
    phenoSimSet = pssObject['data']
    
    # Make sure the model matches in the phenotype simulation set.
    if phenoSimSet['model'] != args.model or phenoSimSet['model_workspace'] != args.workspace:
        print 'Specified model %s/%s does not match model %s/%s in phenotype simulation set' \
            %(args.workspace, args.model, phenoSimSet['model_workspace'], phenoSimSet['model'])
    print "  %d simulations in phenotype simulation set" %(len(phenoSimSet['phenotypeSimulations']))

    # Go through the list of simulations, for each gene and see if the gene is in the genesToReactions
    # dictionary.  If so, see if the gene was a lethal or non-lethal knockout and add it to the 
    # corresponding dictionary along with reactions and 
    lethalList = list()
    nonlethalList = list()
    for sim in phenoSimSet['phenotypeSimulations']:
        if sim[3] == 'CP' or sim[3] == 'CN':
            right = 1
        else:
            right = 0
        geneList = sim[0][0]
        if sim[1] > 0:
            for gene in geneList:
                if gene in genesToReactions:
                    for reaction in genesToReactions[gene]:
                        nonlethalList.append( (gene, reaction, genesToReactions[gene][reaction], right ) )
        else:
            for gene in geneList:
                if gene in genesToReactions:
                    for reaction in genesToReactions[gene]:
                        lethalList.append( (gene, reaction, genesToReactions[gene][reaction], right ) )

    # Sort the lists by reaction probability.
    lethalList.sort(key=itemgetter(2))
    nonlethalList.sort(key=itemgetter(2))
    
    # Walk through each list, generating points for a plot.
    print "  %d knockouts in lethal list" %(len(lethalList))
    right = 0
    wrong = 0
    lethalFile = open(args.phenosimset+'.lethal.csv', 'w')
    lethalFile.write('wrong,right\n')
    for index in range(len(lethalList)):
        if lethalList[index][3]:
            right += 1
        else:
            wrong += 1
        lethalFile.write('%d,%d\n' %(wrong,right))
    lethalFile.close()
    
    print "  %d knockouts in non-lethal list" %(len(nonlethalList))
    right = 0
    wrong = 0
    nonlethalFile = open(args.phenosimset+'.nonlethal.csv', 'w')
    nonlethalFile.write('wrong,right\n')
    for index in range(len(nonlethalList)):
        if nonlethalList[index][3]:
            right += 1
        else:
            wrong += 1
        nonlethalFile.write('%d,%d\n' %(wrong,right))
    nonlethalFile.close()
    exit(0)
    