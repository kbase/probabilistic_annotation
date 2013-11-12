#!/usr/bin/python

# Generate summary statistics comparing two gapfill results tables from AnalyzeGapfillResults.py (one for probanno
# and one for non-probanno)

import optparse
import sys

usage = "%prog [Probanno_result_table] [Non_probanno_result_table]"
description = """Generate summary statistics like number of uniquely-added reactions, etc."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-o", "--added_only", help="Set this flag to ONLY include added reactions (not reversibility changes) in counts of genes (though not reactions) and in average probability calculations.",
                  action="store_true", dest="addedonly", default=False)
parser.add_option("-f", "--final_only", help="Set this flag to only include reactions that are in the final models (e.g. those that were integrated into the model and not deleted by reaction sensitivity analysis) in the counts and average calcualtions. By default we include all reactions in a gapfill solution regardless of downstream analysis", action="store_true", dest="finalonly", default=False)
parser.add_option("-d", "--deleted_only", help="The opposite of -f, set this flag to ONLY print information about gapfilled reactions that have been deleted from the model.", action="store_true", dest="deletedonly", default=False)
parser.add_option("-g", "--print_genes", help="Print detailed gene information instead of summaries.", action="store_true", dest="printgenes", default=False)
(options, args) = parser.parse_args()

if len(args) < 2:
    print usage
    exit(1)

def parse_result_table(filename, addedReactionsOnly=False, finalReactionsOnly=False, deletedOnly = False):
    ''' Parse a results table into a useful data structure keyed by solution number then other information

    '''
    if finalReactionsOnly and deletedOnly:
        raise IOError("Setting both final reactions only and deleted reactions only is a contradiction. Please choose one or the other.")

    fid = open(filename)
    results = {}
    for line in fid:
        # Skip header row
       if "rxn_likelihood" in line:
           continue
       spl = line.strip("\r\n").split("\t")
       solnum = spl[0]
       rxnid = spl[1]
       objval = spl[2]
       gapfill_uuid = spl[3] #not used
       rxn_likelihood = spl[4]
       is_revchange = spl[5]
       gpr = spl[6]
       newgenes = spl[7].split(";")
       numnew = spl[8] #not used (only the aggregate number after taking unique is useful)
       rxn_in_final_model = spl[9]
       if is_revchange == "1" and addedReactionsOnly:
           continue
       if finalReactionsOnly and rxn_in_final_model == "False":
           continue
       if deletedOnly and rxn_in_final_model == "True":
           continue

       if rxn_likelihood == "NO_PROBABILITY":
           rxn_likelihood = "0"
       
       rxninfo = { "likelihood" : rxn_likelihood,
                   "gpr"        : gpr,
                   "newgenes"   : newgenes,
                   "revchange"  : is_revchange }
       if solnum in results:
           results[solnum]["rxninfo"][rxnid] = rxninfo
       else:
           results[solnum] = {}
           results[solnum]["objective"] = objval
           results[solnum]["rxninfo"] = {}
           results[solnum]["rxninfo"][rxnid] = rxninfo
    return results

def getProbabilities(results, solnum, reaction_list, addedReactionsOnly=False):
    ''' Get the probabilities from a set of gapfilled reactions
    '''
    probabilities = []
    for reaction in reaction_list:
        if addedReactionsOnly and results[solnum]["rxninfo"][reaction]["revchange"] == "1":
            continue
        probabilities.append(float(results[solnum]["rxninfo"][reaction]["likelihood"]))
    return probabilities

def getUniqueGenes(results, solnum, reaction_list, addedReactionsOnly = False):
    ''' Get the genes uniquely added from gapfill (e.g. those not in the model previously).
    '''
    unique_genes = set()
    for reaction in reaction_list:
        if addedReactionsOnly and results[solnum]["rxninfo"][reaction]["revchange"] == "1":
            continue
        for gene in results[solnum]["rxninfo"][reaction]["newgenes"]:
            if gene == '':
                continue
            unique_genes.add(gene)
    return unique_genes

def safeAverage(numarray):
    try:
        avg = sum(numarray)/len(numarray)
        return avg
    except ZeroDivisionError:
        return None

probanno_results = parse_result_table(args[0], addedReactionsOnly = options.addedonly, finalReactionsOnly = options.finalonly, deletedOnly = options.deletedonly )
non_probanno_results = parse_result_table(args[1], addedReactionsOnly = options.addedonly, finalReactionsOnly = options.finalonly, deletedOnly = options.deletedonly )

if options.printgenes:
    print "\t".join( [ "reaction", "likelihood", "whenfound", "solnum" ] )
else:
    print "\t".join( [ "probanno_filename", "non_probanno_filename", "solution_number_compared", 
                   "number_common", "number_probanno_only", "number_nonprobanno_only",
                   "average_common", "average_probanno_only", "average_nonprobanno_only",
                   "unique_genes_common", "unique_genes_probanno_only", "unique_genes_nonprobanno_only"]
                 )

for sol in probanno_results.keys():
    if sol not in non_probanno_results:
        continue
    
    # Get reactions in common and unique to each solution
    probanno_reactions     = set(probanno_results[sol]["rxninfo"].keys())
    non_probanno_reactions = set(non_probanno_results[sol]["rxninfo"].keys())
    all_reactions          = probanno_reactions | non_probanno_reactions
    common_reactions       = probanno_reactions & non_probanno_reactions
    unique_to_probanno      = probanno_reactions - non_probanno_reactions
    unique_to_non_probanno = non_probanno_reactions - probanno_reactions

    if options.printgenes:
        for reaction in common_reactions:
            if options.addedonly and probanno_results[sol]["rxninfo"][reaction]["revchange"] == "1":
                continue
            print "%s\t%s\t%s\t%s" %(reaction, str(probanno_results[sol]["rxninfo"][reaction]["likelihood"]), "COMMON", sol)
        for reaction in unique_to_probanno:
            if options.addedonly and probanno_results[sol]["rxninfo"][reaction]["revchange"] == "1":
                continue
            print "%s\t%s\t%s\t%s" %(reaction, str(probanno_results[sol]["rxninfo"][reaction]["likelihood"]), "PROBANNO_ONLY", sol)
        for reaction in unique_to_non_probanno:
            if options.addedonly and non_probanno_results[sol]["rxninfo"][reaction]["revchange"] == "1":
                continue
            print "%s\t%s\t%s\t%s" %(reaction, str(non_probanno_results[sol]["rxninfo"][reaction]["likelihood"]), "NON_PROBANNO_ONLY", sol)
        continue

    # Get unique genes for interesting sets
    common_newgenes                 = getUniqueGenes(probanno_results, sol, common_reactions, addedReactionsOnly = options.addedonly)
    unique_to_probanno_newgenes     = getUniqueGenes(probanno_results, sol, unique_to_probanno, addedReactionsOnly = options.addedonly) - common_newgenes
    unique_to_non_probanno_newgenes = getUniqueGenes(non_probanno_results, sol, unique_to_non_probanno, addedReactionsOnly = options.addedonly) - common_newgenes

    n_common_newgenes = len(common_newgenes)
    n_unique_to_probanno_newgenes = len(unique_to_probanno_newgenes)
    n_unique_to_non_probanno_newgenes = len(unique_to_non_probanno_newgenes)

    # Get probability distributions for interesting sets
    common_probabilities                  = getProbabilities(probanno_results, sol, common_reactions, addedReactionsOnly = options.addedonly)
    unique_to_probanno_probabilities      = getProbabilities(probanno_results, sol, unique_to_probanno, addedReactionsOnly = options.addedonly)
    unique_to_non_probanno_probabilities  = getProbabilities(non_probanno_results, sol, unique_to_non_probanno, addedReactionsOnly = options.addedonly)

    # Get average probabilities for these sets
    common_avg              =  safeAverage(common_probabilities)
    unq_to_probanno_avg     =  safeAverage(unique_to_probanno_probabilities)
    unq_to_non_probanno_avg =  safeAverage(unique_to_non_probanno_probabilities)

    n_common  = len(common_reactions)
    n_unq_to_probanno = len(unique_to_probanno)
    n_unq_to_non_probanno = len(unique_to_non_probanno)

    # Generate historgram data
    # TODO

    print "\t".join( [ args[0], args[1], sol, 
                     str(n_common), str(n_unq_to_probanno), str(n_unq_to_non_probanno), 
                     str(common_avg), str(unq_to_probanno_avg), str(unq_to_non_probanno_avg), 
                     str(n_common_newgenes), str(n_unique_to_probanno_newgenes), str(n_unique_to_non_probanno_newgenes)
                    ] )
