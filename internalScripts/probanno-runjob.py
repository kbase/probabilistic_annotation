#! /usr/bin/python

import sys, os, json
from ConfigParser import ConfigParser
from biokbase.probabilistic_annotation.Impl import ProbabilisticAnnotation

def getConfig(filename):
    retconfig = {}
    config = ConfigParser()
    config.read(filename)
    for nameval in config.items("ProbabilisticAnnotation"):
        retconfig[nameval[0]] = nameval[1]
    return retconfig

# First parameter is the path to the file containing the job object json.
# Read the job object from the file.
job = json.load(open(sys.argv[1], "r"))

# Second parameter is the path to the config file.
# Read the config from the file.
config = getConfig(sys.argv[2])

# Call the server function to run the annotation.
impl_ProbabilisticAnnotation = ProbabilisticAnnotation(config)
impl_ProbabilisticAnnotation.runAnnotate(job["jobdata"])

# Temporary code to clean up.
os.remove(sys.argv[1])

exit(0)
