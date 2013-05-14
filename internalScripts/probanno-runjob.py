#! /usr/bin/python

import sys, os, json
from biokbase.probabilistic_annotation.Impl import ProbabilisticAnnotation
from biokbase.probabilistic_annotation.DataParser import getConfig

# First parameter is the path to the file containing the job object json.
# Read the job object from the file.
job = json.load(open(sys.argv[1], "r"))

# Second parameter is the path to the config file.
# Read the config from the file.
config = getConfig(sys.argv[2])

# Let the server know that it is being called from a job.
config["generate_data_option"] = "runjob"

# Call the server function to run the annotation.
impl_ProbabilisticAnnotation = ProbabilisticAnnotation(config)
impl_ProbabilisticAnnotation.runAnnotate(job["jobdata"])

# Temporary code to clean up.
os.remove(sys.argv[1])

exit(0)
