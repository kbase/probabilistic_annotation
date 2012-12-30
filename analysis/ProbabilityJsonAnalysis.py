#!/usr/bin/python
try:
    import json
except ImportError:
    sys.path.append('simplejson-2.3.3')
    import simplejson as json

import os
import operator
import sys
import re

class probabilityJson:
    def __init__(self, f):
        if not os.path.exists(f):
            raise IOError
        self.json = json.load(open(f, "r"))
        self.genedict = self.makeGeneIdDict()
        # TC is for transporters.
        self.ecfinder = re.compile("\(\s*[ET]\.*C\.*\s*[\-0-9]+\.[\-0-9]+\.[\-0-9]+\.[\-0-9]+\s*\)")
    def improveAnnotationConsistency(self, function):
        '''Attempt to improve consistency between different sets of roles
        by fixing two very common problems:
        1: EC number in one role but not the other, and
        2: Capitalization issues.
        '''
        noec = self.ecfinder.sub("", function)
        lowercase = noec.lower()
        # Remove extraneous white space
        stripped = lowercase.strip()
        return stripped
    def makeGeneIdDict(self):
        '''Make a dictionary from gene ID to feature object'''
        features=self.json["features"]
        id2feature = {}
        for feature in features:
            if feature["type"] == "peg" or feature["type"] == "CDS":
                id2feature[feature["id"]] = feature
        return id2feature
    def listPegIds(self):
        return self.genedict.keys()
    def getGeneFunction(self, fid):
        '''Get the function (from KBase) for the specified gene'''
        if fid not in self.genedict:
            raise KeyError("Specified gene ID %s not found in the JSON file" %(fid))
        feature = self.genedict[fid]
        if "function" not in feature:
            return None
        elif feature["function"] == "":
            return None
        else:
            return self.improveAnnotationConsistency(feature["function"])
    def getProbabilisticFunctions(self, fid):
        '''For a single feature ID, create a dictionary from gene function predicted by
        probabilitstic annotation to the probability of that function'''
        if fid not in self.genedict:
            raise KeyError("Specified gene ID %s not found in the JSON file" %(fid))
        feature = self.genedict[fid]
        # No probabilities could be calculated based on the source data set
        if "alternativeFunctions" not in feature:
            return None
        else:
            func2prob={}
            for af in feature["alternativeFunctions"]:
                improvedFunction = self.improveAnnotationConsistency(af[0])
                func2prob[improvedFunction] = af[1]
            return func2prob

s = probabilityJson(sys.argv[1])
ls = s.listPegIds()

# Assessment of agreement between probabilistic and
# the functional annotation in the JSON file...
probNoReal = 0
realNoProb = 0
realNoProbHypothetical = 0
realNoProbAmbiguous = 0
noEither = 0
bothAgree = 0
doNotAgree = 0
agreeAnnoteToProb = {}
maxDisagreeAnnoteToProb = {}
disagreeRealToOthers = {}
disagreeRealToMaxAnnote = {}
for fid in ls:
    probfunc = s.getProbabilisticFunctions(fid)
    realfunc = s.getGeneFunction(fid)
    # No function either from the genome object OR from probabilities.
    if probfunc is None and realfunc is None:
        noEither += 1
        continue
    if probfunc is None and realfunc is not None:
        realNoProb += 1
        if "hypothetical" in realfunc or "Hypothetical" in realfunc or "unknown" in realfunc or "Uncharacterized" in realfunc or "uncharacterized" in realfunc:
            realNoProbHypothetical += 1
        elif "putative" in realfunc or "Putative" in realfunc:
            realNoProbAmbiguous += 1
        continue
    if realfunc is None and probfunc is not None:
        probNoReal += 1
        continue
    if realfunc in probfunc:
        bothAgree += 1
        agreeAnnoteToProb[realfunc] = probfunc[realfunc]
    else:
#        print realfunc
#        print probfunc
#        print ""

        doNotAgree += 1
        mxkey = max(probfunc.iteritems(), key=operator.itemgetter(1))
        maxDisagreeAnnoteToProb[mxkey[0]] = mxkey[1]
        disagreeRealToOthers[realfunc] = "\t".join( [ k + "(" + str(probfunc[k]) + ")" for k in probfunc ] )
        disagreeRealToMaxAnnote[realfunc] = mxkey[0] + "\t" + str(mxkey[1])

print "NoEither: %d" %(noEither)
print "realNoProb: %d" %(realNoProb)
print "realNoProbHypothetical: %d" %(realNoProbHypothetical)
print "realNoProbAmbiguous: %d" %(realNoProbAmbiguous)
print "probNoReal: %d" %(probNoReal)
print "bothAgree: %d" %(bothAgree)
print "doNotAgree: %d" %(doNotAgree)

#for func in disagreeRealToOthers:
#    print "%s\t%s" %(func, disagreeRealToOthers[func])

for func in disagreeRealToMaxAnnote:
    print "%s\t%s" %(func, disagreeRealToMaxAnnote[func])

#print "agreeAnnoteToProb:"
#print agreeAnnoteToProb
#print "maxDisagreeAnnoteToProb:"
#print maxDisagreeAnnoteToProb
