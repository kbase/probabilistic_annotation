import os, sys, tempfile, shutil
from biokbase.probabilistic_annotation.DataExtractor import *
from biokbase.probabilistic_annotation.DataParser import *
from biokbase.probabilistic_annotation.PYTHON_GLOBALS import *
from biokbase.workspaceService.client import *
from biokbase.fbaModelServices.Client import *
import scipy
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import numpy

def genomeToFasta(params, genomeObject, workFolder):
    '''Convert a Genome object into an amino-acid FASTA file (for BLAST purposes)'''
    
    # Open the fasta file.
    fastaFile = os.path.join(workFolder, "%s.faa" %(params.genome))
    try:
        fid = open(fastaFile, "w")
    except IOError:
        # TODO Throw an exception back to the client
        print("Failed to open fasta file %s\n" %(fastaFile))
        return
    
    # Make sure the Genome object has features.
    # TODO Throw exception if no features in genome
    if "features" not in genomeObject["data"]:
         sys.stderr.write("""ERROR: The input Genome object %s/%s has no features - did you forget to run annotate_genome?\n""" %(params.genome_workspace, params.genome))
         return None

    # Run the list of features to build the fasta file.
    features = genomeObject["data"]["features"]
    for feature in features:
        # Not a protein-encoding gene
        if "protein_translation" not in feature:
            continue
        myid = feature["id"]
        if "function" in feature:
            function = feature["function"]
        else:
            function = ""
        seq = feature["protein_translation"]
        fid.write(">%s %s\n%s\n" %(myid, function, seq))
        
    return fastaFile
    
def runBlast(params, queryFile, workFolder, dataFolder):
    '''A simplistic wrapper to BLAST the query proteins against the subsystem proteins'''
    
    blastResultFile = os.path.join(workFolder, "%s.blastout" %(params.genome))
    cmd = "blastp -query \"%s\" -db %s -outfmt 6 -evalue 1E-5 -num_threads 1 -out \"%s\"" %(queryFile, os.path.join(dataFolder, SUBSYSTEM_OTU_FASTA_FILE), blastResultFile)
    sys.stderr.write("Started BLAST with command: %s\n" %(cmd))
    status = os.system(cmd)
    sys.stderr.write("Ended BLAST with command: %s\n" %(cmd))
    # TODO Throw exception if blast fails
    if status != 0:
        sys.stderr.write("BLAST command failed")
    return blastResultFile

def rolesetProbabilitiesMarble(params, blastResultFile, workFolder, dataFolder):
    '''Calculate the probabilities of rolesets (i.e. each possible combination of
    roles implied by the functions of the proteins in subsystems) from the BLAST results.

    Returns a file with three columns(query, rolestring, probability)  
    rolestring = "\\\" separating all roles of a protein with a single function (order does not matter)
    '''

    sys.stderr.write("Performing marble-picking on rolesets...")

    # Read in the target roles (this function returns the roles as lists!)
    targetIdToRole, targetRoleToId = readFilteredOtuRoles(dataFolder)

    # Convert the lists of roles into "rolestrings" (sort the list so that order doesn't matter)
    # in order to deal with the case where some of the hits are multi-functional and others only have
    # a single function...
    targetIdToRoleString = {}
    for target in targetIdToRole:
        stri = SEPARATOR.join(sorted(targetIdToRole[target]))
        targetIdToRoleString[target] = stri

    # Query --> [ (target1, score 1), (target 2, score 2), ... ]
    idToTargetList = parseBlastOutput(blastResultFile)

    # This is a holder for all of our results
    # It is a list of tuples in the form (query, rolestring, probability)
    rolestringTuples = {}
    # For each query gene we calcualte the likelihood of each possible rolestring.
    for query in idToTargetList:
        # First we need to know the maximum score
        # I have no idea why but I'm pretty sure Python is silently turning the second element of these tuples
        # into strings.
        #
        # That's why I turn them back...
        maxscore = 0
        for tup in idToTargetList[query]:
            if float(tup[1]) > maxscore:
                maxscore = float(tup[1])

        # Now we calculate the cumulative squared scores
        # for each possible rolestring. This along with PC*maxscore is equivalent
        # to multiplying all scores by themselves and then dividing by the max score.
        # This is done to avoid some pathological cases and give more weight to higher-scoring hits
        # and not let much lower-scoring hits \ noise drown them out.
        rolestringToScore = {}
        for tup in idToTargetList[query]:
            try:
                rolestring = targetIdToRoleString[tup[0]]
            except KeyError:
#                sys.stderr.write("ERROR: Target ID %s from the BLAST file had no roles in the rolestring dictionary??\n" %(tup[0]))
                continue
            if rolestring in rolestringToScore:
                rolestringToScore[rolestring] += (float(tup[1]) ** 2)
            else:
                rolestringToScore[rolestring] = (float(tup[1])**2)

        # Now lets iterate over all of them and calculate the probability
        # Probability = sum(S(X)^2)/(sum(S(X)^2 + PC*maxscore))
        # Lets get the denominator first
        denom = PSEUDOCOUNT*maxscore
        for stri in rolestringToScore:
            denom += rolestringToScore[stri]

        # Now the numerators, which are different for every rolestring
        for stri in rolestringToScore:
            p = rolestringToScore[stri]/denom
            if query in rolestringTuples:
                rolestringTuples[query].append( (stri, p) )
            else:
                rolestringTuples[query] = [ (stri, p) ]

    # Save the generated data when debug is turned on.
    if params.debug:
        rolesetProbabilityFile = os.path.join(workFolder, "%s.rolesetprobs" %(params.genome))
        fid = open(rolesetProbabilityFile, "w")
#        for f, p in rolestringTuples.iteritems():
        for query in rolestringTuples:
            for tup in rolestringTuples[query]:
                fid.write("%s\t%s\t%1.4f\n" %(query, tup[0], tup[1]))
        fid.close()
        
    sys.stderr.write("done\n")
    return rolestringTuples
        
#def MakeProbabilisticJsonFile(annotation_file, blast_result_file, roleset_probability_file, outfile, folder, genome_id, probanno_id):
def buildProbAnnoObject(params, genomeObject, blastResultFile, queryToRolesetProbs, workFolder, dataFolder, wsClient):
    '''Create a "probabilistic annotation" object file from a Genome object file. The probabilistic annotation
    object adds fields for the probability of each role being linked to each gene.'''

    sys.stderr.write("Building ProbAnno object for genome %s..." %(params.genome))
    targetToRoles, rolesToTargets = readFilteredOtuRoles(dataFolder)
    targetToRoleSet = {}
    for target in targetToRoles:
        stri = SEPARATOR.join(sorted(targetToRoles[target]))
        targetToRoleSet[target] = stri
    queryToTargetProbs = parseBlastOutput(blastResultFile)
    
    # For each query ID:
    # 1. Identify their rolestring probabilities (these are the first and second elements of the tuple)
    # 2. Iterate over the target genes and identify those with each function (a list of these and their blast scores is
    #    the third element of the tuple) - should be able to set it up to only iterate once over the list.
    # 3. Add that tuple to the JSON file with the key "alternativeFunctions"

    # The Genome object data ["features"] is a list of dictionaries. We want to make our data structure and 
    # then add that to the dictionary.  I use the ii in range so I can edit the elements without changes being lost.

    objectData = {}
    objectData["id"] = params.probanno
    objectData["genome"] = params.genome
    objectData["genome_workspace"] = params.genome_workspace;
    objectData["rolesetProbabilities"] = queryToRolesetProbs;
    objectData["featureAlternativeFunctions"] = []
    objectData["skippedFeatures"] = []
    
    for ii in range(len(genomeObject["data"]["features"])):
        feature = genomeObject["data"]["features"][ii]
        if "id" not in genomeObject["data"]:
            # TODO Throw an exception
            sys.stderr.write("No gene ID found in input Genome object (this should never happen)\n")
            return None
        queryid = feature["id"]

        # This can happen if I couldn't find hits from that gene to anything in the database. In this case, I'll just skip it.
        # TODO Or should I make an empty object? I should ask Chris.
        if queryid not in queryToRolesetProbs or queryid not in queryToTargetProbs:
            objectData["skippedFeatures"].append(queryid)
            continue
        
        # Get a list of (target ID, BLAST score) pairs for each possible role set from the query's homologs
        # At the moment we dont' use this...
        # I'm keeping the code around in case we decide we want the functionality.
        query_rolesetToTargetScores = {}
        for targetprob in queryToTargetProbs[queryid]:
            try:
                targetroleset = targetToRoleSet[targetprob[0]]
            except KeyError:
                sys.stderr.write("WARNING: Internal data inconsistency - BLAST target %s has no associated roleset\n" %(targetprob[0]))
                continue
            if targetroleset in query_rolesetToTargetScores:
                query_rolesetToTargetScores[targetroleset].append(targetprob)
            else:
                query_rolesetToTargetScores[targetroleset] = [ targetprob ]

        # For each possible roleset, obtain its probability, and the list of (target, BLAST score) pairs
        featureAlternativeFunctions = {}
        alternative_functions = []
        for rolesetprob in queryToRolesetProbs[queryid]:
            roleset = rolesetprob[0]
            rolep = rolesetprob[1]
            if roleset not in query_rolesetToTargetScores:
                # TODO Throw exception
                sys.stderr.write("INTERNAL ERROR: The rolesets for the query genes were not transferred from targets correctly - query roleset %s did not exist in the target homolog list\n" %(roleset))
                exit(2)
            rlist = [ roleset, rolep ]
            alternative_functions.append(rlist)
        featureAlternativeFunctions["alternativeFunctions"] = alternative_functions
        featureAlternativeFunctions["id"] = queryid
        objectData["featureAlternativeFunctions"].append(featureAlternativeFunctions)

    # Store the ProbAnno object in the specified workspace.
    objectMetadata = { "num_rolesets": len(objectData["rolesetProbabilities"]),
                       "num_altfuncs": len(objectData["featureAlternativeFunctions"]),
                       "num_skipped_features": len(objectData["skippedFeatures"]) }
    saveObjectParams = { "type": "ProbAnno", "id": params.probanno, "workspace": params.probanno_workspace,
                         "auth": params.auth, "data": objectData, "metadata": objectMetadata, "command": "pa-annotate" }
    metadata = wsClient.save_object(saveObjectParams)
    
    sys.stderr.write("done\n")
    return metadata

def rolesetProbabilitiesToRoleProbabilities(params, queryToTuplist, workFolder):
    '''Compute probability of each role from the rolesets for each query protein.
    At the moment the strategy is to take any set of rolestrings containing the same roles
    and add their probabilities.
    So if we have hits to both a bifunctional enzyme with R1 and R2, and
    hits to a monofunctional enzyme with only R1, R1 ends up with a greater
    probability than R2.

    I had tried to normalize to the previous sum but I need to be more careful than that
    (I'll put it on my TODO list) because if you have e.g.
    one hit to R1R2 and one hit to R3 then the probability of R1 and R2 will be unfairly
    brought down due to the normalization scheme...

    Returns an array where each element has three columns: Query gene ID, role, and probability'''

    sys.stderr.write("Generating role probabilities from roleset probabilities...")

    roleProbs = []    
    for query in queryToTuplist:
        # TODO MATT It looks like totalp is not used in this version of the code
        # Get the total probability - needed for re-normalization later.
#         totalp = 0
#         for tup in queryToTuplist[query]:
#             totalp += tup[1]

        # This section actually does the conversion of probabilities.
        queryRolesToProbs = {}
        for tup in queryToTuplist[query]:
            rolelist = tup[0].split(SEPARATOR)
            # Add up all the instances of each particular role on the list.
            for role in rolelist:
                if role in queryRolesToProbs:
                    queryRolesToProbs[role] += tup[1]
                else:
                    queryRolesToProbs[role] = tup[1]

        # Add them to the array.
        for role in queryRolesToProbs:
            roleProbs.append( (query, role, queryRolesToProbs[role]) )

    # Save the generated data when debug is turned on.
    if params.debug:
        role_probability_file = os.path.join(workFolder, "%s.roleprobs" %(params.genome))
        fid = open(role_probability_file, "w")
        for tuple in roleProbs:
            fid.write("%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2]))
        fid.close()

    sys.stderr.write("done\n")
    return roleProbs

# For now to get the probability I just assign this as the MAXIMUM for each role
# to avoid diluting out by noise.
#
# The gene assignments are all genes within DILUTION_PERCENT of the maximum...
#
def totalRoleProbabilities(params, roleProbs, workFolder):
    '''Given the probability that each gene has each role, estimate the probability that
    the entire ORGANISM has that role.

    To avoid exploding the probabilities with noise, I just take the maximum probability
    of any query gene having a function and use that as the probability that the function
    exists in the cell.

    Returns a file with three columns: each role, its probability, and the estimated set of genes
    that perform that role. A gene is assigned to a role if it is within DILUTION_PERCENT
    of the maximum probability. DILUTION_PERCENT is defined in the python_globals.py file.
    '''

    sys.stderr.write("Generating whole-cell role probability file...")
    

    # Compute maximum probability among all query genes for each role.
    roleToTotalProb = {}
    for tuple in roleProbs:
        if tuple[1] in roleToTotalProb:
            if float(tuple[2]) > roleToTotalProb[tuple[1]]:
                roleToTotalProb[tuple[1]] = float(tuple[2])
        else:
            roleToTotalProb[tuple[1]] = float(tuple[2])

    # Get the genes within DILUTION_PERCENT percent of the maximum
    # probability and assert those (note - DILUTION_PERCENT is defined in the global parameters file)
    # This is a dictionary from role to a list of genes
    roleToGeneList = {}
    for tuple in roleProbs:
        if tuple[1] not in roleToTotalProb:
            # TODO Throw an exception
            sys.stderr.write("ERROR: Role %s not placed properly in roleToTotalProb dictionary?\n" %(tuple[1]))
            return None
        if float(tuple[2]) >= float(DILUTION_PERCENT)/100.0 * roleToTotalProb[tuple[1]]:
            if tuple[1] in roleToGeneList:
                roleToGeneList[tuple[1]].append(tuple[0])
            else:
                roleToGeneList[tuple[1]] = [ tuple[0] ]
    
    # Build the array of total role probabilities.     
    totalRoleProbs = []
    for role in roleToTotalProb:
        totalRoleProbs.append( (role, roleToTotalProb[role], " or ".join(roleToGeneList[role])) )   

    # Save the generated data when debug is turned on.
    if params.debug:
        total_role_probability_file = os.path.join(workFolder, "%s.cellroleprob" %(params.genome))
        fid = open(total_role_probability_file, "w")
        for tuple in totalRoleProbs:
            fid.write("%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2]))
        fid.close()
        
    sys.stderr.write("done\n")
    return totalRoleProbs

def complexProbabilities(params, totalRoleProbs, workFolder, dataFolder):
    '''Compute the probability of each complex from the probability of each role.

    The complex probability is computed as the minimum probability of roles within that complex.
    
    The output file to this has four columns
    Complex   |   Probability   | Type   |  Roles_not_in_organism  |  Roles_not_in_subsystems
    
    Type: 

    CPLX_FULL (all roles found and utilized)
    CPLX_PARTIAL (only some roles found - only those roles that were found were utilized; does not distinguish between not there and no reps for those not found)
    CPLX_NOTTHERE (Probability of 0 because the genes aren't there for any of the subunits)
    CPLX_NOREPS (Probability 0f 0 because there are no representative genes in the subsystems)
    '''

    sys.stderr.write("Computing complex probabilities...")

    # Get the mapping from complexes to roles.
    complexesToRequiredRoles = readComplexRoles(dataFolder)
    
    # Get the subsystem roles (used to distinguish between NOTTHERE and NOREPS).
    otu_fidsToRoles, otu_rolesToFids  = readFilteredOtuRoles(dataFolder)
    allroles = set()
    for fid in otu_fidsToRoles:
        for role in otu_fidsToRoles[fid]:
            allroles.add(role)

    # 2: Read the total role --> probability file
    rolesToProbabilities = {}
    rolesToGeneList = {}
    for tuple in totalRoleProbs:
        rolesToProbabilities[tuple[0]] = float(tuple[1]) # can skip the float()?
        rolesToGeneList[tuple[0]] = tuple[2]

    # Iterate over complexes and compute complex probabilities from role probabilities.
    # Separate out cases where no genes seem to exist in the organism for the reaction from cases
    # where there is a database deficiency.
    complexProbs = []
    for cplx in complexesToRequiredRoles:
        allCplxRoles = complexesToRequiredRoles[cplx]
        availRoles = [] # Roles that may have representatives in the query organism
        unavailRoles = [] # Roles that have representatives but that are not apparently in the query organism
        noexistRoles = [] # Roles with no representatives in the subsystems
        for role in complexesToRequiredRoles[cplx]:
            if role not in allroles:
                noexistRoles.append(role)
            elif role not in rolesToProbabilities:
                unavailRoles.append(role)
            else:
                availRoles.append(role)
        TYPE = ""
        GPR = ""
        if len(noexistRoles) == len(allCplxRoles):
            TYPE = "CPLX_NOREPS"
            complexProbs.append( (cplx, 0.0, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR) )
#            fid.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(cplx, 0.0, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR))
            continue
        if len(unavailRoles) == len(allCplxRoles):
            TYPE = "CPLX_NOTTHERE"
            complexProbs.append( (cplx, 0.0, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR) )
#            fid.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(cplx, 0.0, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR))
            continue
        # Some had no representatives and the rest were not found in the cell
        if len(unavailRoles) + len(noexistRoles) == len(allCplxRoles):
            TYPE = "CPLX_NOREPS_AND_NOTTHERE"
            complexProbs.append( (cplx, 0.0, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR) )
#            fid.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(cplx, 0.0, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR))
            continue
        # Otherwise at least one of them is available
        if len(availRoles) == len(allCplxRoles):
            TYPE = "CPLX_FULL"
        elif len(availRoles) < len(allCplxRoles):
            TYPE = "CPLX_PARTIAL_%d_of_%d" %(len(availRoles), len(allCplxRoles))
        GPR = " and ".join("(" + s + ")" for s in [ rolesToGeneList[f] for f in availRoles ] )
        if GPR != "==":
            GPR = "(" + GPR + ")"
        # Find the minimum probability of the different available roles (ignoring ones that are apparently missing)
        # and call that the complex probability
        minp = 1000
        for role in availRoles:
            if rolesToProbabilities[role] < minp:
                minp = rolesToProbabilities[role]
        complexProbs.append( (cplx, minp, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR) )
#        fid.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(cplx, minp, TYPE, SEPARATOR.join(unavailRoles), SEPARATOR.join(noexistRoles), GPR))

    if params.debug:
        complex_probability_file = os.path.join(workFolder, "%s.complexprob" %(params.genome))
        fid = open(complex_probability_file, "w")
        for tuple in complexProbs:
            fid.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4], tuple[5]))
        fid.close()
    
    sys.stderr.write("done\n")
    return complexProbs

def reactionProbabilities(params, complexProbs, workFolder, dataFolder):
    '''From the probability of complexes estimate the probability of reactions.

    The reaction probability is computed as the maximum probability of complexes that perform
    that reaction.

    The output to this function is a file containing the following columns:
    Reaction   |  Probability  |  RxnType  |  ComplexInfo

    If the reaction has no complexes it won't even be in this file becuase of the way
    I set up the call... I could probably change this so that I get a list of ALL reactions
    and make it easier to catch issues with reaction --> complex links in the database.
    Some of the infrastructure is already there (with the TYPE).

    ComplexInfo is information about the complex IDs, their probabilities, and their TYPE
    (see ComplexProbabilities)
    '''

    sys.stderr.write("Computing reaction probabilities...")
    
    # Cplx --> {Probability, type, GPR}
    cplxToTuple = {}
    for tuple in complexProbs:
        cplxToTuple[tuple[0]] = ( tuple[1], tuple[2], tuple[5] )
    
    # Take the MAXIMUM probability of complexes
    rxnsToComplexes = readReactionComplex(dataFolder)

    reactionProbs = []
    for rxn in rxnsToComplexes:
        TYPE = "NOCOMPLEXES"
        rxnComplexes = rxnsToComplexes[rxn]
        complexinfo = ""
        maxp = 0
        GPR = ""
        for cplx in rxnComplexes:
            if cplx in cplxToTuple:
                # Complex1 (P1; TYPE1) ///Complex2 (P2; TYPE2) ...
                complexinfo = "%s%s%s (%1.4f; %s) " %(complexinfo, SEPARATOR, cplx, cplxToTuple[cplx][0], cplxToTuple[cplx][1])
                TYPE = "HASCOMPLEXES"
                if cplxToTuple[cplx][0] > maxp:
                    maxp = cplxToTuple[cplx][0]
                if GPR == "":
                    GPR = cplxToTuple[cplx][2]
                elif cplxToTuple[cplx][2] != "":
                    GPR = " or ".join( [ GPR, cplxToTuple[cplx][2] ] )
        reactionProbs.append( (rxn, maxp, TYPE, complexinfo, GPR) )

    if params.debug:
        reaction_probability_file = os.path.join(workFolder, "%s.rxnprobs" %(params.genome))
        fid = open(reaction_probability_file, "w")
        for tuple in reactionProbs:
            fid.write("%s\t%1.4f\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]))
        fid.close()

    # (end of function)
    sys.stderr.write("done\n")
    return reactionProbs

def buildModelObject(params, reactionProbs):
    ''' Create a probability model object from the reaction probabilities.'''
    
    sys.stderr.write("Building model object...")
    
    # Build the list of reaction probabilities.
    rxnProbList = []
    for rxn in reactionProbs:
        rxnProbList.append( (rxn[0], rxn[1], rxn[4]) )
        
    # Create a fba modeling service client.
    fbaModelClient = fbaModelServices(FBAMODEL_URL)
    
    # Build a model from the input genome and reaction probabilities.
    probFbaModelParams = { "genome": params.genome, "genome_workspace": params.genome_workspace,
                           "model": params.model, "workspace": params.model_workspace,
                           "reaction_probs": rxnProbList, "default_prob": 0.0,
                           "overwrite": params.overwrite, "auth": params.auth }
    metadata = fbaModelClient.genome_to_probfbamodel(probFbaModelParams)
    
    sys.stderr.write("done\n")
    return metadata

def metaboliteWeights(params, model):
    
    '''
    Given a model object with name "model",
    computes an S-matrix.
 
    This function returns three things:
    - A sparse matrix object (coo_matrix) from scipy with indexes matching the lists
    - A dictionary from metabolite UUIDs to their index in the matrix
    - A dictionary from reaction UUIDs to their index in the matrix
 
    The lists should have the same order as the input model.
 
    If absval = True, returns the absolute value of the S-matrix (absolute value of every
    term in the S matrix) rather than S itself.
     
    TODO - put in biomass equation too.
    '''
 
    print(params)
    metIdToIdx = {}
    rxnIdToIdx = {}
 
    idx = 0
    for compound in model["modelcompounds"]:
        metIdToIdx[compound["uuid"]] = idx
        idx += 1
     
    i = []
    j = []
    data = []
    costs = []
    idx = 0
    numZeroReagents = 0
    for reaction in model["modelreactions"]:
        if "modelReactionReagents" not in reaction:
            numZeroReagents += 1
            continue
        for reagent in reaction["modelReactionReagents"]:
            met_uuid = reagent["modelcompound_uuid"]
            coefficient = reagent["coefficient"]
            met_idx = metIdToIdx[met_uuid]
            i.append(met_idx)
            j.append(idx)
            if params.absval:
                coefficient = abs(coefficient)
            data.append(coefficient)
        costs.append(1.0 - reaction["probability"])
        rxnIdToIdx[reaction["uuid"]] = idx
        idx += 1
 
    S_matrix = sparse.coo_matrix( ( data, ( i, j ) ) )
    
#     for biomass in model["biomasses"]:
#         for cpd in biomass["biomasscompounds"]:
#             sys.stderr.write("%s: %d\n" %(cpd["modelcompound_uuid"], metIdToIdx[cpd["modelcompound_uuid"]]))
                  
    sys.stderr.write("%d model reactions had zero reagents\n" %(numZeroReagents))
    sys.stderr.write("len i=%d len j=%d len data=%d\n" %(len(i), len(j), len(data)))
    sys.stderr.write("S_matrix rows=%d, cols=%d\n" %(S_matrix.shape[0], S_matrix.shape[1]))

    rxncosts = numpy.array(costs)    
    '''
    Given an S matrix (sparse) S, and a vector of reaction 
    probabilities rxnprobs, calls the scipy least-squares solver
    to obtain our best estimate for the metabolite weights
    (gamma)
 
    The reaction weights (probabilities) and metabolite weights
    are related by the equation
    R = |S|^T * gamma
 
    where R is the vector of reaction weights, |S| means the absolute value
    of S and gamma is the vector of metabolite weights. Solving in the least-squares
    sense is the best we can do since S is not a square matrix.
 
    Returns a list of metabolite weights with the same indexing as the rows of S.
    '''
     
    S_prime = S_matrix.transpose(copy=True)
    res = linalg.lsqr(S_prime, rxncosts)
    
#    for index in range(len(model["modelcompounds"])):
    for index in range(len(res[0])):
        cpd = model["modelcompounds"][index]
        cpd["weight"] = res[0][index]
        
    return True

def safeRemove(fname, dirname):
    totalfname = os.path.join(dirname, fname)
    try:
        # Check for file existence
        fid = open(totalfname, "r")
        fid.close()
        os.remove(totalfname)
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass

