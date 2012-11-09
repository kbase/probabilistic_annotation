#!/usr/bin/python

import os
import tempfile

URL="http://bio-data-1.mcs.anl.gov/services/cdmi_api"

# This is a character string that isn't found in any roles - I use it to split lists
# in one-to-many relationships
SEPARATOR = "///"

# NOTE - This will be the location of the data folder
# from the workspace - CHRIS: Edit this when we have
# the workspace working! I have no idea where you plan
# on putting this directory when you make it...
BASEDIR=tempfile.gettempdir()

if not os.path.isdir(os.path.join(BASEDIR, "data")):
    os.path.mkdir(os.path.join(BASEDIR, "data"))

# Names of temporary files containing all of the data
# that is NOT query-specific
SUBSYSTEM_FID_FILE = os.path.join(BASEDIR, "data", "SUBSYSTEM_FIDS")
DLIT_FID_FILE = os.path.join(BASEDIR, "data", "DLIT_FIDS")
CONCATINATED_FID_FILE = os.path.join(BASEDIR, "data", "ALL_FIDS")
OTU_ID_FILE = os.path.join(BASEDIR, "data", "OTU_GENOME_IDS")
SUBSYSTEM_OTU_FIDS_FILE = os.path.join(BASEDIR, "data", "SUBSYSTEM_OTU_FIDS")
SUBSYSTEM_OTU_FID_ROLES_FILE = os.path.join(BASEDIR, "data", "SUBSYSTEM_OTU_FID_ROLES")
SUBSYSTEM_OTU_FASTA_FILE = os.path.join(BASEDIR, "data", "SUBSYSTEM_FASTA")
#OTU_NEIGHBORHOOD_FILE = "data/OTU_NEIGHBORHOODS"
COMPLEXES_ROLES_FILE = os.path.join(BASEDIR, "data", "COMPLEXES_ROLES")
REACTION_COMPLEXES_FILE = os.path.join(BASEDIR, "data", "REACTIONS_COMPLEXES")

# Minimum index
MINN = 0
# Number of items to grab
COUNT = 50000

# Percentage of the maximum probability to use as a threshold
# to consider other genes as having a particular function aside
# from the one with greatest probability
#
# D: 80% (will play with this)
DILUTION_PERCENT = 80
