#!/usr/bin/python

import argparse
import sys
import re
import os
import traceback
from biokbase.probabilistic_annotation.Client import _read_inifile
from biokbase.workspaceService.Client import *
from biokbase.fbaModelServices.Client import *
from biokbase.cdmi.client import CDMI_API, CDMI_EntityAPI

desc = '''
Given a JSON object calculate the false positive, negative rates. Optionally omit everything that isn't
in the specified media.
'''
separator ='///'
rolesToRxnsFile = 'roles_to_reactions'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='AnalyzePhenotypeSimulationResults.py', description=desc)
    parser.add_argument('genome', help='ID of genome to analyze', action='store', default=None)
    parser.add_argument('model', help='ID of integrated gap filled model to analyze', action='store', default=None)
    parser.add_argument('phenosimset', help='ID of PhenotypeSimulationSet object to analyze', action='store', default=None)
    parser.add_argument('workspace', help='ID of workspace containing objects to analyze', action='store', default=None)
    parser.add_argument('--solution', help='index number of solution in integrated gap filled model', action='store', default='0')
    parser.add_argument('--phenosimsetws', help='workspace containing PhenotypeSimulationSet object to analyze', action='store', default=None)
    parser.add_argument('--media', help='limit analysis to only this media condition', action='store', dest='media', default=None)
    parser.add_argument('--probanno', help='ID of ProbAnno object with probabilities', action='store', dest='probanno', default=None)
    parser.add_argument('--probannows', help='workspace containing ProbAnno object (same as phenosimsetws if not specified)', action='store', dest='probannows', default=None)
    parser.add_argument('--script-dir', help='path to directory with analysis scripts', action='store', dest='scriptDir', default='.')
    parser.add_argument('--fba-url', help='url for fba modeling service', action='store', dest='fbaurl', default='http://bio-data-1.mcs.anl.gov/services/fba')
    parser.add_argument('--ws-url', help='url for workspace service', action='store', dest='wsurl', default='http://www.kbase.us/services/workspace/')
    args = parser.parse_args()
    
    if args.probanno is not None and args.probannows is None:
        args.probannows = args.phenosimsetws
    if args.phenosimsetws is None:
        args.phenosimsetws = args.workspace
        
    # Create the clients for the services we need.
    authdata = _read_inifile()
    token = authdata['token']
    wsClient = workspaceService(args.wsurl)
    
    # Build roles to reactions dictionary.
    step = 1
    print "+++ Step %d: Build reactions to roles dictionary +++" %(step)
    if not os.path.exists(rolesToRxnsFile):
        os.system("all_roles_used_in_models | roles_to_complexes | get_relationship_HasStep -to id >%s" %(rolesToRxnsFile))
    rolesToReactions = dict()
    reactionsToRoles = dict()
    # Each line of the file has four fields: (1) role, (2) hypothetical flag, (3) complex id, (4) reaction id
    handle = open(rolesToRxnsFile, 'r')
    for line in handle:
        fields = line.strip('\r\n').split('\t')
        reaction = fields[3]
        if reaction not in reactionsToRoles:
            reactionsToRoles[reaction] = list()
        reactionsToRoles[reaction].append(fields[0])
#        rolesToReactions[fields[0]] = { 'hypothetical': fields[1], 'complex': fields[2], 'reaction': fields[3] }
    handle.close()
    print "  %d reactions in reactions to roles dictionary" %(len(reactionsToRoles))

    # Analyze the gap fill results for the specified model.
    step += 1
    print "+++ Step %d: Run AnalyzeGapfillResults.py for model '%s/%s' +++" %(step, args.workspace, args.model)
    rxnprobs = args.genome + '.rxnprobs'
    resultsFile = args.model + '.results'
    os.system("%s/AnalyzeGapfillResults.py -m %s -w %s --rxnprobs %s --url %s >%s" \
              %(args.scriptDir, args.model, args.workspace, rxnprobs, args.fbaurl, resultsFile))
#    reactionDict = dict()
    genesToReactions = dict()
    handle = open(resultsFile, 'r')
    handle.readline() # Throw away header line
    for line in handle:
        fields = line.strip('\r\n').split('\t')
        if fields[0] == args.solution:
            if fields[5] == '0': # Only keep reactions that are not a reversibility change
                if fields[6]: # Only keep reactions that have a GPR
                    rxnid = re.sub(r'rxn0*(\d+)', r'kb|rxn.\1', fields[1])
                    geneList = re.findall('fig\|\d+\.\d+\.peg\.\d+', fields[6])
                    for index in range(len(geneList)):
#                        print rxnid + ': ' + geneList[index]
                        if geneList[index] not in genesToReactions:
                            genesToReactions[geneList[index]] = dict()
                        genesToReactions[geneList[index]][rxnid] = 0.0
#                    reactionDict[rxnid] = { 'likelihood': fields[4], 'gpr': fields[6], 'newgenes': fields[7], 'numnewgenes': fields[8] }
    handle.close()
    print "  %d genes in genes to reactions dictionary" %(len(genesToReactions))
#    print "  %d reactions in reaction dictionary" %(len(reactionDict))

    # Get the ProbAnno object.
    step += 1
    probanno = args.genome + '.probanno'
    print "+++ Step %d: Get ProbAnno object '%s/%s'" %(step, args.workspace, probanno)
    paObject = wsClient.get_object( { 'id': probanno, 'workspace': args.workspace, 'type': 'ProbAnno', 'auth': token } )
    probAnno = paObject['data']
    print "  %d genes in ProbAnno roleset probabilities dictionary" %(len(probAnno['roleset_probabilities']))
    
    # Need to go through rolesets dictionary
    # If an entry has more than one role, split into parts
    # Then run the array and combine any duplicate roles by adding probs
    step += 1
    print "+++ Step %d: Parse rolesets into roles and adjust probabilities for duplicates +++" %(step)
    rolesetProbabilities = dict()
    for gene in probAnno['roleset_probabilities']:
        geneRoleList = probAnno['roleset_probabilities'][gene]
        geneRoleDict = dict()
        for index in range(len(geneRoleList)):
            prob = geneRoleList[index][1]
            roleList = geneRoleList[index][0].split(separator)
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
    probsFile = open(args.genome+'probs', 'w')
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
    print "+++ Step %d: Get phenotype simulation set %s/%s +++" %(step, args.phenosimset, args.phenosimsetws)
    pssObject = wsClient.get_object( { 'id': args.phenosimset, 'workspace': args.phenosimsetws, 'type': 'PhenotypeSimulationSet', 'auth': token } )
    phenoSimSet = pssObject['data']
    
    # Make sure the model matches in the phenotype simulation set.
    if phenoSimSet['model'] != args.model or phenoSimSet['model_workspace'] != args.workspace:
        print 'Specified model %s/%s does not match model %s/%s in phenotype simulation set' \
            %(args.workspace, args.model, phenoSimSet['model_workspace'], phenoSimSet['model'])
    print "  %d simulations in phenotype simulation set" %(len(phenoSimSet['phenotypeSimulations']))

    step += 1
    print '+++ Step %d: Calculate curve +++' %(step) 
    csvFile = open(args.phenosimset+'.csv', 'w')
    csvFile.write('Media,Cutoff,TN,TP,FN,FP\n')
    prob = 0.0
    while prob <= 1.0:
        mediaDict = dict()
        for sim in phenoSimSet['phenotypeSimulations']:
            media = sim[0][1]
            if media not in mediaDict:
                mediaDict[media] = dict()
                mediaDict[media]['trueNeg'] = 0
                mediaDict[media]['truePos'] = 0
                mediaDict[media]['falseNeg'] = 0
                mediaDict[media]['falsePos'] = 0
            geneList = sim[0][0]
            if len(geneList) > 1:
                print 'got more than 1 gene in genelist'
            keep = False
            for gene in geneList:
                if gene in genesToReactions:
                    for reaction in genesToReactions[gene]:
                        if genesToReactions[gene][reaction] >= prob:
                            keep = True
#                else:
#                    print 'gene %s in sim set is not in geneToReactions dict' %(gene)
            if keep:
                simType = sim[3]
                if simType == "CN":
                    mediaDict[media]['trueNeg'] += 1
                elif simType == "CP":
                    mediaDict[media]['truePos'] += 1
                elif simType == "FN":
                    mediaDict[media]['falseNeg'] += 1
                elif simType == "FP":
                    mediaDict[media]['falsePos'] += 1
        prob += 0.01
        for media in mediaDict:
            csvFile.write("%s,%f,%d,%d,%d,%d\n" %(media, prob, mediaDict[media]['trueNeg'], mediaDict[media]['truePos'], mediaDict[media]['falseNeg'], mediaDict[media]['falsePos']))
    csvFile.close()
    exit(0)
    