#!/usr/bin/python

# Main wrapper script for front-end calculations

import optparse
import os, sys
from DataProcessor import *

usage="%prog [options] organismid"
description="""Main driver to run probabilistic annotation.

First generates a list of Query-independent (i.e. same for all queries)
data from the KBase ER model. This includes lists of OTUs (pre-computed
as part of the KBase) and the members of subsystems that are in those OTUs;
also included are their roles.

Then, it uses the provided organism ID to look for a JSON file and, if
it does not exist already, generates one for you from the central store (note
that the function to do this for now actually generates a new base genome ID for you
- this is ignored)...

Finally, it does a BLAST against the pre-computed data, and uses that result
to come up with annotation probabilities.
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--regenerate", help="Regenerate database if it already exists (NOTE - takes a long time)", 
                  action="store_true", dest="regenerate", default=False)
(options, args) = parser.parse_args()

# Make folders where results will go...
try:
    os.mkdir("data")
except OSError:
    pass;

try:
    os.mkdir(os.path.join("data", "OTU"))
except OSError:
    pass;

if len(args) < 1:
    sys.stderr.write("ERROR: Organism ID is a required argument\n")
    exit(2)

# Run the extractor driver to get the data
cmd = "python ExtractorDriver.py"
if options.regenerate:
    cmd += cmd + " -r"
os.system(cmd)

organismid = args[0]
fasta_file, json_file = setUpQueryData(organismid)
blast_result_file = runBlast(organismid, fasta_file)
roleset_probability_file = RolesetProbabilitiesMarble(organismid, blast_result_file)
role_probability_file = RolesetProbabilitiesToRoleProbabilities(organismid, roleset_probability_file)
total_role_probability_file = TotalRoleProbabilities(organismid, role_probability_file)
complex_probability_file = ComplexProbabilities(organismid, total_role_probability_file)
reaction_probability_file = ReactionProbabilities(organismid, complex_probability_file)

outfile = os.path.join(organismid, "%s_prob.json" %(organismid))
MakeProbabilisticJsonFile(json_file, blast_result_file, roleset_probability_file, outfile)