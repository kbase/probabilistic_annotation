#!/usr/bin/python

import optparse, sys
from biokbase.probabilistic_annotation.MyImpl import generate_data

usage="%prog [options]"
description="""Main driver to get data needed out of the KBase and store it locally.
Data will be stored in a local database autoReconInfo
All of this data is QUERY INDEPENDENT. It should all be the same
for any organism for which you want to do a reconstruction..."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--regenerate", help="Regenerate database if it already exists (NOTE - takes a long time)", action="store_true", dest="regenerate", default=False)
parser.add_option("-d", "--deleteonly", help="Delete data files but do not regenerate (WARNING - this is not reversible)", action="store_true", dest="delete_only", default=False)
parser.add_option("-v", "--verbose", help="Display all WARNINGS (D: Only display messages related to completeness)", action="store_true", dest="verbose", default=False)
parser.add_option("-f", "--folder", help="Base directory (folder) in which all of the data files are to be stored", action="store", dest="folder", default=None)
(options, args) = parser.parse_args()

if options.folder is None:
    sys.stderr.write("ERROR: In ExtractorDriver.py - folder (-f) is a required argument\n")
    exit(2)

success = generate_data(options)
if success:
    exit(0)
exit(1)
