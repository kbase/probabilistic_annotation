#!/usr/bin/python

# This program will apply pseudocounts to calculate a probability that a protein has each of several specified functions.
# The "marble picking scheme" implemented here does the following:
#
# Take the log-likelihood of each of the E-values
# Parse them out according to what function (NOT ROLE) they belong to 
#
# Calculate the probability of a gene having a particular function as:
# [Sum for FUNCTION X/Total elements of FUNCTION X] / [ [Sum_X]/[N_FUNCTIONX] + [Sum_Y]/[N_FUNCTIONY] + ... + PSEUDOCOUNT]
#
# The PSEUDOCOUNT is set to avoid: 0/0
# and also to avoid one very weak hit getting 100% probability
# It is set so that if the there is exactly one hit with an E-value of 1E-10 (meaning the sum of means is 10)
# then the probability of that function is 10/(10+PC) = 50% or PC = 10
#
# The means are to normalize for the number of things in the database (since the functions are quite divergent in terms of
# E-value we can't just normalize by number of things with a given function)
#
######
#
# Input looks like this:
#
# BLAST_results: standard -outfmt 6 result of a blastp run
#
# TargetFunctions: Two-column table
# gene   |  function
#
# The list of genes must be the same as the list of genes in the BLAST database.
#
##### 
# Output comes out like this:
#
# gene | function | probability
#
# These probabilities are NOT normalized by number

import sys
import re
import math
import functools

E_PSEUDOCOUNT=10
N_STDEV = 3
#S_PERCENTILE=0.7
MIN_EVALUE=1E-200
MAX_EVALUE=1
PROB_CUTOFF = 0.01 # Don't bother reporting probabilities lower than this...

if not len(sys.argv) == 3:
    print "Usage: ./BlastResToPseudocounts.py [BLAST_results] [TargetFunctions]"
    exit(3)

class TargetPair(object):
    def __init__(self, targetid, evalue):
        self.targetid = str(targetid)
        self.evalue = float(evalue)
    def getEvalue(self):
        return self.evalue
    def getTarget(self):
        return self.targetid
    def getPair(self):
        return self.targetid, self.evalue

# Container for lists of (Target, Evalue) pairs
# "Target" could be a FIGFAM, a fucntion, whatever.
class TargetPairList(object):
    def __init__(self, TP = None):
        if TP is None:
            self.TPList = []
        else:
            self.TPList = [ TP ]
    def addTargetPair(self, targetid, evalue):
        self.TPList.append(TargetPair(targetid, evalue))
    def getTargetList(self):
        # Get a list of targets in the current list.
        TargetList = []
        for TP in self.TPList:
            TargetList.append(TP.getTarget())
        return list(set(TargetList))
    def getEvalueList(self):
        # Get a list of all of the E-values in the current list.
        EvalueList = []
        for TP in self.TPList:
            EvalueList.append(TP.getEvalue())
        return EvalueList
    def getEvaluesForTarget(self, targetid):
        # Returns an empty array if no entries with the target are found.
        EvalueList = []
        for TP in self.TPList:
            if TP.getTarget() == targetid:
                EvalueList.append(TP.getEvalue())
        return EvalueList

def getNormalizedLog(E, mine, maxe):
    if E < mine:
        E = mine
    if E > maxe:
        E = maxe
    return -1.0 * math.log10(E)

# Treat stdev of a lengh-1 array as 0 for purpose of our calculation
def getStdev(arr):
    if len(arr) == 1:
        return 0
    mean = sum(arr)/len(arr)
    sumsq = sum( [ (i-mean)**2 for i in arr ] )
    return (sumsq / (len(arr) - 1) )**0.5


def percentile(N, percent, key=lambda x:x):
    """
    From http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

# median is 50th percentile.
median = functools.partial(percentile, percent=0.5)

# Mapping from (target) gene ID to the functions (not to be confused with functional roles)
geneToFunction = {}
sys.stderr.write("Reading Target Functions...\n")
for line in open(sys.argv[2], "r"):
    spl = line.strip().split("\t")
    geneToFunction[spl[0]] = spl[1]
    
# Mapping from a query gene ID to a list of target \ E-value pairs.
queryToTPList = {}

sys.stderr.write("Reading BLAST results...\n")
for line in open(sys.argv[1], "r"):
    spl = line.strip().split("\t")
    # Initialize a new set of target-Evalue pairs if necessary
    if not spl[0] in queryToTPList:
        queryToTPList[spl[0]] = TargetPairList()
    if not spl[1] in geneToFunction:
        sys.stderr.write("WARNING: No function found for gene %s\n" %(spl[1]))
        queryToTPList[spl[0]].addTargetPair("UNKNOWN", spl[10])
    else:
        queryToTPList[spl[0]].addTargetPair(geneToFunction[spl[1]] , spl[10])

# Finally, we iterate over the queries, calculate the sum of log-likelihoods for
# each of the targets in it and divide by the appropriate normalization factors.
sys.stderr.write("Calculating probabilities...\n")
for query in queryToTPList:
    TargetList = queryToTPList[query].getTargetList()
    EvalueList = queryToTPList[query].getEvalueList()
    means = []

    # Calculate the sum of the means [+ Pseudocount to avoid cases where 1 high hit to a figfam would wash out lots of high hits to a different figfam]
    for target in TargetList:
        TargetEvalues = queryToTPList[query].getEvaluesForTarget(target)
        LogTargetEvalues = []
        for E in TargetEvalues:
            LogTargetEvalues.append(getNormalizedLog(E, MIN_EVALUE, MAX_EVALUE))
        # Mean + k Stdev
        mean = sum(LogTargetEvalues)/len(LogTargetEvalues)
        means.append(sum(LogTargetEvalues)/len(LogTargetEvalues) + N_STDEV*getStdev(LogTargetEvalues) )
        #means.append(mean)
        #means.append(percentile(sorted(LogTargetEvalues), S_PERCENTILE))

    meansum = sum(means)

    for target in TargetList:
        TargetEvalues = queryToTPList[query].getEvaluesForTarget(target)
        LogTargetEvalues = []
        for E in TargetEvalues:
            LogTargetEvalues.append(getNormalizedLog(E, MIN_EVALUE, MAX_EVALUE))
        # Add up the mean + k * stdev
        prob = ( sum(LogTargetEvalues)/len(LogTargetEvalues) + N_STDEV*getStdev(LogTargetEvalues) )/(meansum + E_PSEUDOCOUNT)
        #prob = sum(LogTargetEvalues)/len(LogTargetEvalues) / (meansum + E_PSEUDOCOUNT)
        #prob = percentile(sorted(LogTargetEvalues), S_PERCENTILE)/(meansum + E_PSEUDOCOUNT)

        if prob < PROB_CUTOFF:
            continue

        print "%s\t%s\t%s" %(query, target, prob)
