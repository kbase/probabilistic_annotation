#!/usr/bin/python

import json
import sys

'''

This is very much a thrown-together script. Sometime i'll set down and do a better job integrating this
into the modelseed API.

Given a JSON object calculate the false positive, negative rates. Optionally omit everything that isn't
in the specified media.

'''

if len(sys.argv) < 1:
    raise IOError("Usage: $prog [Phenotype_json_file] [Media]")

s = json.load(open(sys.argv[1], "r"))

if len(sys.argv) > 2:
    media = sys.argv[2]
else:
    media = None

sims = s["phenotypeSimulations"]

CN = 0
CP = 0
FN = 0
FP = 0

for sim in sims:
    med = sim[0][1]
    typ = sim[3]
    if media is not None:
        if media != med:
            continue
    if typ == "CN":
        CN +=1
    elif typ == "CP":
        CP += 1
    elif typ == "FN":
        FN += 1
    elif typ == "FP":
        FP += 1

print "CN (correct negative): %d" %(CN)
print "CP (correct positive): %d" %(CP)
print "FN (false negative): %d" %(FN)
print "FP (false positive): %d" %(FP)
