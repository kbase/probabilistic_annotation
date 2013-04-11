#!/usr/bin/python

import os
import tempfile

CDMI_URL="http://kbase.us/services/cdmi_api"

# This is a character string that isn't found in any roles - I use it to split lists
# in one-to-many relationships
SEPARATOR = "///"

# Names of temporary files containing all of the data
# that is NOT query-specific
SUBSYSTEM_FID_FILE = "SUBSYSTEM_FIDS"
DLIT_FID_FILE = "DLIT_FIDS"
CONCATINATED_FID_FILE = "ALL_FIDS"
CONCATINATED_FID_ROLE_FILE = "ALL_FID_ROLES"
OTU_ID_FILE = "OTU_GENOME_IDS"
SUBSYSTEM_OTU_FID_ROLES_FILE = "SUBSYSTEM_OTU_FID_ROLES"
SUBSYSTEM_OTU_FASTA_FILE = "SUBSYSTEM_FASTA"
#OTU_NEIGHBORHOOD_FILE = "data/OTU_NEIGHBORHOODS"
COMPLEXES_ROLES_FILE = "COMPLEXES_ROLES"
REACTION_COMPLEXES_FILE = "REACTIONS_COMPLEXES"

# Minimum index
MINN = 0
# Number of items to grab
#COUNT = 50000
COUNT = 50000

# Percentage of the maximum probability to use as a threshold
# to consider other genes as having a particular function aside
# from the one with greatest probability
#
# D: 80% (will play with this)
DILUTION_PERCENT = 80
