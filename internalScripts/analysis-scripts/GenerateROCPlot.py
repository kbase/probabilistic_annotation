#!/usr/bin/python -u

import argparse
import sys
import re
import os
import traceback
from operator import itemgetter, attrgetter
import matplotlib.pyplot as plt
import numpy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='GenerateROCPlot.py', description='make a plot')
    parser.add_argument('fileList', help='list of results file names', nargs='+', action='store', default=None)
    parser.add_argument('--outputFile', help='file name', action='store', default='output.pdf')
    args = parser.parse_args()

    totalPos = 0
    totalNeg = 0
    resultList = list()
    for name in args.fileList:
        resultsFile = open(name, "r")
        resultsFile.readline() # Throw away header line
        for line in resultsFile:
            fields = line.strip('\n').split(',')
            resultList.append( [ float(fields[2]), fields[3] ] )
            if fields[3] == 'CP' or fields[3] == 'FP':
                totalPos += 1
            else:
                totalNeg += 1
    
    print "totalPos %d totalNeg %d" %(totalPos, totalNeg)
    resultList.sort(key=itemgetter(0), reverse=True)
    rocFile = open(args.outputFile, 'w')
    rocFile.write("prob,totalPos,numTP,numFP,tpr,fpr,totalNeg,numTN,numFN,tnr,fnr\n")
    for cutoff in numpy.arange(0,1,0.05):
        numTP = 0
        numFP = 0
        numTN = 0
        numFN = 0
        for i in range(len(resultList)):
            if resultList[i][0] >= cutoff:
                if resultList[i][1] == 'CP':
                    numTP += 1
                elif resultList[i][1] == 'FP':
                    numFP += 1
                elif resultList[i][1] == 'TN':
                    numTN += 1
                else:
                    numFN += 1
        if totalPos > 0:
            tpr = float(numTP)/float(totalPos)
            fpr = float(numFP)/float(totalPos)
        else:
            tpr = float(0.0)
            fpr = float(0.0)
        if totalNeg > 0:
            tnr = float(numTN)/float(totalNeg)
            fnr = float(numFN)/float(totalNeg)
        else:
            tnr = float(0.0)
            fnr = float(0.0)
        rocFile.write("%f,%d,%d,%d,%f,%f,%d,%d,%d,%f,%f\n" %(cutoff, totalPos, numTP, numFP, tpr, fpr, totalNeg, numTN, numFN, tnr, fnr) )
    rocFile.close()

#     plt.plot(tpr, fpr)
#     plt.savefig(args.outputFile)
