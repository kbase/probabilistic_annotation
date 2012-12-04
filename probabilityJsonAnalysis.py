
try:
    import json
except:
    import simplejson as json

import os
import operator
import sys

class probabilityJson:
    def __init__(self, f):
        if not os.path.exists(f):
            raise IOError
        self.json = json.load(open(f, "r"))
        self.genedict = self.makeGeneIdDict()
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
            return feature["function"]
    def getProbabilisticFunctions(self, fid):
        '''For a single feature ID, create a dictionary from gene function to 
        the probabilistic annotations'''
        if fid not in self.genedict:
            raise KeyError("Specified gene ID %s not found in the JSON file" %(fid))
        feature = self.genedict[fid]
        # No probabilities could be calculated based on the source data set
        if "alternativeFunctions" not in feature:
            return None
        else:
            func2annote={}
            for af in feature["alternativeFunctions"]:
                func2annote[af[0]] = af[1]
            return func2annote

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
        doNotAgree += 1
        mxkey = max(probfunc.iteritems(), key=operator.itemgetter(1))
        maxDisagreeAnnoteToProb[mxkey[0]] = mxkey[1]
        disagreeRealToOthers[realfunc] = "\t".join( [ k + "(" + str(probfunc[k]) + ")" for k in probfunc ] )

print "NoEither: %d" %(noEither)
print "realNoProb: %d" %(realNoProb)
print "realNoProbHypothetical: %d" %(realNoProbHypothetical)
print "realNoProbAmbiguous: %d" %(realNoProbAmbiguous)
print "probNoReal: %d" %(probNoReal)
print "bothAgree: %d" %(bothAgree)
print "doNotAgree: %d" %(doNotAgree)

for func in disagreeRealToOthers:
    print "%s\t%s" %(func, disagreeRealToOthers[func])
#print "agreeAnnoteToProb:"
#print agreeAnnoteToProb
#print "maxDisagreeAnnoteToProb:"
#print maxDisagreeAnnoteToProb
