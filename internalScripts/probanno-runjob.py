#! /usr/bin/python

import sys, os, json
from biokbase.probabilistic_annotation.Impl import ProbabilisticAnnotation
from biokbase.probabilistic_annotation.DataParser import getConfig

# First parameter is the path to the job directory.
# Read the job object from the file in the job directory.
jsonFilename = os.path.join(sys.argv[1], "jobfile.json")
job = json.load(open(jsonFilename, "r"))

# Read the config from the file specified in the job data.
config = job["jobdata"]["config"]

# Let the server know that it is being called from a job.
config["load_data_option"] = "runjob"

# Call the server function to run the annotation.
impl_ProbabilisticAnnotation = ProbabilisticAnnotation(config)
impl_ProbabilisticAnnotation.runAnnotate(job)

exit(0)
