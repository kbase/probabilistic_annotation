#!/usr/bin/python

import json
import argparse
import sys
import traceback
from biokbase.probabilistic_annotation.Client import _read_inifile
from biokbase.fbaModelServices.Client import *
from biokbase.workspace.client import *

desc = '''
Given a JSON object calculate the false positive, negative rates. Optionally omit everything that isn't
in the specified media.
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='AnalyzePhenotypeSimulationResults.py', description=desc)
    parser.add_argument('pheno', help='ID of PhenotypeSimulationSet object to analyze', action='store', default=None)
    parser.add_argument('phenows', help='Workspace containing PhenotypeSimulationSet object to analyze', action='store', default=None)
    parser.add_argument('--media', help='Limit analysis to only this media', action='store', dest='media', default=None)
    parser.add_argument('--rxnprobs', help='ID of RxnProbs object with probabilities for reactions associated with false positive/negative genes', action='store', dest='rxnprobs', default=None)
    parser.add_argument('--rxnprobsws', help='workspace containing RxnProbs object (same as phenows if not specified)', action='store', dest='rxnprobsws', default=None)
    parser.add_argument('--csv', help='Generate output in csv format', action='store_true', dest='csv', default=False)
    parser.add_argument('--split', help='Split results by media', action='store_true', dest='split', default=False)
    parser.add_argument('--fba-url', help='url for fba modeling service', action='store', dest='fbaurl', default='https://kbase.us/services/KBaseFBAModeling')
    parser.add_argument('--ws-url', help='url for workspace service', action='store', dest='wsurl', default='http://kbase.us/services/ws/')
    parser.add_argument('--cdmi-url', help='url for cdmi service', action='store', dest='cdmiurl', default='http://kbase.us/services/cdmi_api/')
    args = parser.parse_args()

    if args.rxnprobs is not None and args.rxnprobsws is None:
        args.rxnprobsws = args.phenows

    # Get simulation data
    wsClient = Workspace(args.wsurl)
    phenoIdentity = { 'workspace': args.phenows, 'name': args.pheno }
    object_list = wsClient.get_objects( [ phenoIdentity ] )
    phenoobj = object_list[0]


    # We need to get reaction IDs associated with each gene.
    model_ref = phenoobj["data"]["fbamodel_ref"]
    model_objectIdentity = { 'ref' : model_ref }
    model_info = wsClient.get_object_info( [ model_objectIdentity ], 0 )
    modelid = model_info[0][1]
    modelws = model_info[0][7]

    fbaClient = fbaModelServices(args.fbaurl)
    models = fbaClient.get_models( { "models" : [ modelid ],
                                     "workspaces" : [ modelws ]
                                     })
    model = models[0]
    modeldict = {}
    for reaction in model["reactions"]:
        features = reaction["features"]
        for feature in features:
            if feature in modeldict:
                modeldict[feature].append( [ reaction["reaction"], reaction["gapfilled"] ] )
            else:
                modeldict[feature] = [ [ reaction["reaction"], reaction["gapfilled"] ] ]

    # Get the RxnProbs object and Model object if requested.
    rxnprobdict = {}

    if args.rxnprobs is not None:
        # Get the RxnProbs object.
        rxnprobsIdentity = { 'workspace' : args.rxnprobsws, 'name': args.rxnprobs }
        object_list = wsClient.get_objects( [ rxnprobsIdentity ] )
        rxnprobobj = object_list[0]
        # RXNPROBS: reaction -> rxnprob_array ( [ rxnid, probability, diagnostic, complexes, GPR ] )
        for rxnprob in rxnprobobj["data"]["reaction_probabilities"]:
            rxnprobdict[rxnprob[0]] = rxnprob
            pass
        pass

    if not args.csv:
        if args.media is None:
            print "For phenotypedata %s and all media: " %(args.pheno)
        else:
            print "For phenotypedata %s and media %s: " %(args.pheno, args.media)

    # Get the phenotype dataset object to which the phenotype simulation refers.
    # We need this to tell what gene(s) were knocked out and what media was used to
    # do the simulation.
    phenotypeIdentity = { 'ref' : phenoobj["data"]["phenotypeset_ref"] }
    object_list = wsClient.get_objects( [ phenotypeIdentity ] )
    phenodata = object_list[0]
    pheno_id_to_info = {}
    for pheno in phenodata["data"]["phenotypes"]:
        pheno_id_to_info[pheno["id"]] = ( pheno["media_ref"], pheno["geneko_refs"] )

    # Count the CP\CN\FP\FN
    # If a rxnprobs object is specified, print the rxnprobs data associated with each gene...
    mediaDict = dict()
    # need to track CP/CN/FP/FN by media
    for sim in phenoobj["data"]["phenotypeSimulations"]:

        # Reach into phenotype data object to get media and gene KO info
        spl = sim["phenotype_ref"].split("/")
        phenoid = spl[-1]
        if phenoid not in pheno_id_to_info:
            raise IOError("Phenotype ID %s not found in the linked phenotype dataset" %(phenoid) )

        media_ref  = pheno_id_to_info[phenoid][0]
        media_ObjectIdentity = { 'ref' : media_ref }
        media_info = wsClient.get_object_info( [ { 'ref' : media_ref } ], 0 ) 
        media = media_info[0][1]

        gene_ko_refs = pheno_id_to_info[phenoid][1]
        genes = []
        for ref in gene_ko_refs:
            spl = ref.split("/")
            genes.append(spl[-1])

        # Some genes are not in the model
        ok = 1
        for gene in genes:
            if gene not in modeldict:
                ok = 0
                break
        if ok == 0:
#            sys.stderr.write("One of these genes not in model: %s\n" %("\t".join(genes)))
            continue

        # Initialize results for a specific media.
        if media not in mediaDict:
            mediaDict[media] = dict()
            mediaDict[media]['trueNeg'] = 0
            mediaDict[media]['truePos'] = 0
            mediaDict[media]['falseNeg'] = 0
            mediaDict[media]['falsePos'] = 0
            mediaDict[media]['unknown'] = 0
        sim_type = sim["phenoclass"]
        if args.media is not None:
            if args.media != media:
                continue

        if len(rxnprobdict.keys()) > 0:
            for gene in genes:
                for rxnpair in modeldict[gene]:
                    # Some reactions in the model have no probabilities associated (we give them 0)
                    rxn = rxnpair[0]
                    isgapfilled = str(rxnpair[1])
                    if rxn in rxnprobdict:
                        rxnprobobj = rxnprobdict[rxn]
                    else:
                        rxnprobobj = [ rxn, 0, "NOCOMPLEX", "", "" ]
                    # Gene, WasGapfilled, reaction, probability, diagnostics...
                    print "\t".join( [ args.pheno, gene, media, isgapfilled ] + [ str(s) for s in rxnprobobj ] + [ sim_type ] )

        if sim_type == "CN":
            mediaDict[media]['trueNeg'] += 1
        elif sim_type == "CP":
            mediaDict[media]['truePos'] += 1
        elif sim_type == "FN":
            mediaDict[media]['falseNeg'] += 1
        elif sim_type == "FP":
            mediaDict[media]['falsePos'] += 1
        elif sim_type == "UN":
            mediaDict[media]['unknown'] += 1
            

    # Calculate sensitivity and specificity for each media.
    for media in mediaDict:
        if (mediaDict[media]['truePos'] + mediaDict[media]['falseNeg']) > 0:
            mediaDict[media]['sensitivity'] = float(mediaDict[media]['truePos']) / float(mediaDict[media]['truePos'] + mediaDict[media]['falseNeg'])
        else:
            mediaDict[media]['sensitivity'] = 0.0
        if (mediaDict[media]['falsePos'] + mediaDict[media]['trueNeg']) > 0:
            mediaDict[media]['specificity'] = float(mediaDict[media]['trueNeg']) / float(mediaDict[media]['falsePos'] + mediaDict[media]['trueNeg'])
        else:
            mediaDict[media]['specificity'] = 0.0
        mediaDict[media]['fpr'] = 1.0 - mediaDict[media]['specificity']

    # Set type based on suffix in phenotype simulation set name.
    if '.std.' in args.pheno:
        type = 'std'
    elif '.pa.' in args.pheno:
        type = 'pa'
    else:
        type = 'unknown'

    # Report results
    if args.split:
        for media in mediaDict:
            if args.csv:
                print "%s,%s,%s,%f,%f,%f,%d,%d,%d,%d" %(args.pheno, type, media, mediaDict[media]['fpr'], \
                                                        mediaDict[media]['sensitivity'], mediaDict[media]['specificity'], \
                                                        mediaDict[media]['trueNeg'], mediaDict[media]['truePos'], \
                                                        mediaDict[media]['falseNeg'], mediaDict[media]['falsePos'])
            else:
                print "Media: %s" %(media)
                print "False positive rate: %f" %(mediaDict[media]['fpr'])
                print "Sensitivity: %f" %(mediaDict[media]['sensitivity'])
                print "Specificity: %f" %(mediaDict[media]['specificity'])
                print "True negative: %d" %(mediaDict[media]['trueNeg'])
                print "True positive: %d" %(mediaDict[media]['truePos'])
                print "False negative: %d" %(mediaDict[media]['falseNeg'])
                print "False positive: %d" %(mediaDict[media]['falsePos'])
    else:
        trueNeg = 0
        truePos = 0
        falseNeg = 0
        falsePos = 0
        for media in mediaDict:
            trueNeg += mediaDict[media]['trueNeg']
            truePos += mediaDict[media]['truePos']
            falseNeg += mediaDict[media]['falseNeg']
            falsePos += mediaDict[media]['falsePos']
        if (truePos + falseNeg) > 0:
            sensitivity = float(truePos) / float(truePos + falseNeg)
        else:
            sensitivity = 0.0
        if (falsePos + trueNeg) > 0:
            specificity = float(trueNeg) / float(falsePos + trueNeg)
        else:
            specificity = 0.0
        fpr = 1.0 - specificity
        if args.csv:
            print "%s,%s,%s,%f,%f,%f,%d,%d,%d,%d" %(args.pheno, type, media, fpr, sensitivity, specificity, \
                                                    trueNeg, truePos, falseNeg, falsePos)
        else:
            print "False positive rate: %f" %(fpr)
            print "Sensitivity: %f" %(sensitivity)
            print "Specificity: %f" %(specificity)
            print "True negative: %d" %(trueNeg)
            print "True positive: %d" %(truePos)
            print "False negative: %d" %(falseNeg)
            print "False positive: %d" %(falsePos)

    exit(0)
