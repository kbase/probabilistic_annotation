#!/usr/bin/python

'''
These are functions for parsing a ModelSEED Model object
and for converting reaction probabilities into metabolite
probabilities.

The conversion is done by using the least-squares approach
( minimize ||b-Ax||^2 ) using scipy's linalg.lsqr() function
'''

import scipy
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

def S_from_model(model, absval = False):
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

    metIdToIdx = {}
    rxnIdToIdx = {}

    idx = 0
    for compound in model["modelcompounds"]:
        metIdToIdx[compound["uuid"]] = idx
        idx += 1
    
    i = []
    j = []
    data = []
    idx = 0
    for reaction in model["modelreactions"]:
        for reagent in reaction["modelReactionReagents"]:
            met_uuid = reagent["modelcompound_uuid"]
            coefficient = reagent["coefficient"]
            met_idx = metIdToIdx[met_uuid]
            i.append(met_idx)
            j.append(idx)
            if absval:
                coefficient = abs(coefficient)
            data.append(coefficient)
        rxnIdToIdx[reaction["uuid"]] = idx
        idx += 1


    matrix = sparse.coo_matrix( ( data, ( i, j ) ) )
    return matrix, metIdToIdx, rxnIdToIdx

def rxnProbsFromModel(model):
    '''
    Get a list of reaction probabilities from a model object.

    The list will be in the same order as the input.
    '''

    return probs

def rxnProbsToMetProbs(S, rxnprobs):
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
    
    S_prime = S.transpose(copy=True)
    res = linalg.lsqr(S_prime, rxnprobs)
    # res[0] is the actual computed value of gamma...
    return res[0]

import json
import sys
model = json.load(open(sys.argv[1], "r"))
res = S_from_model(model, absval=False)
rxnprobs = rxnprobs_from_model(model)
print res[0]
