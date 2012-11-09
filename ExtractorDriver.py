#!/usr/bin/python

import optparse, os, sqlite3, sys
from DataExtractor import *
from DataParser import *
# Common variables...
from PYTHON_GLOBALS import *

usage="%prog [options]"
description="""Main driver to get data needed out of the KBase and store it locally.
Data will be stored in a local database autoReconInfo
All of this data is QUERY INDEPENDENT. It should all be the same
for any organism for which you want to do a reconstruction..."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--regenerate", help="Regenerate database if it already exists (NOTE - takes a long time)", action="store_true", dest="regenerate", default=False)
parser.add_option("-d", "--deleteonly", help="Delete data files but do not regenerate (WARNING - this is not reversible)", action="store_true", dest="delete", default=False)
parser.add_option("-v", "--verbose", help="Display all WARNINGS (D: Only display messages related to completeness)", action="store_true", dest="verbose", default=False)
(options, args) = parser.parse_args()

def safeRemove(f):
    try:
        # Check for file existence
        fid = open(f, "r")
        fid.close()
        os.remove(f)
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass

if options.regenerate or options.delete:
    safeRemove(OTU_ID_FILE)
    safeRemove(SUBSYSTEM_FID_FILE)
    safeRemove(DLIT_FID_FILE)
    safeRemove(CONCATINATED_FID_FILE)
    safeRemove(SUBSYSTEM_OTU_FIDS_FILE)
    safeRemove(SUBSYSTEM_OTU_FID_ROLES_FILE)
    safeRemove(SUBSYSTEM_OTU_FASTA_FILE)
    safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".psq") 
    safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".pin")
    safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".phr")
    safeRemove(OTU_NEIGHBORHOOD_FILE)
    safeRemove(COMPLEXES_ROLES_FILE)
    safeRemove(REACTION_COMPLEXES_FILE)

#    folder = os.path.join("data", "OTU")
#    for the_file in os.listdir(folder):
#        file_path = os.path.join(folder, the_file)
#        if os.path.isfile(file_path):
#            os.unlink(file_path)

# Our job is done if all we want to do is delete files.
if options.delete:
    exit(0)

sys.stderr.write("Generating requested data:....\n")

############
# Get lists of OTUs
############
sys.stderr.write("OTU data...\n")
try:
    otus, prokotus = readOtuData()
except IOError:
    otus, prokotus = getOtuGenomeIds(MINN, COUNT)
#    otus, prokotus = getOtuGenomeIds(MINN, 1200)
    writeOtuData(otus, prokotus)

############
# Get a list of subsystem FIDs
############
sys.stderr.write("List of subsystem FIDS...\n")
try:
    sub_fids = readSubsystemFids()
except IOError:
    sub_fids = subsystemFids(MINN, COUNT)
    # NOTE - This is a TEMPORARY workaround for an issue with
    # the KBase subsystem load. This function WILL BE DELETED
    # and reverted to the call above once that issue is fixed...
#    sub_fids = subsystemFids_WORKAROUND(MINN, COUNT)
    writeSubsystemFids(sub_fids)

###########
# ALso get a list of Dlit FIDs
# We include these because having them
# greatly expands the number of roles for which we
# have representatives.
##########
sys.stderr.write("Getting a list of DLit FIDs...\n")
try:
    dlit_fids = readDlitFids()
except IOError:
    dlit_fids = getDlitFids(MINN, COUNT)
    writeDlitFids(dlit_fids)

##########
# Concatinate the two FID lists before filtering
# (Note - doing so after would be possible as well but
# can lead to the same kinds of biases as not filtering
# the subsystems... Im not sure the problem would
# be as bad for these though)
##########
sys.stderr.write("Combining lists of subsystem and DLit FIDS...\n")
try:
    all_fids = set()
    for line in open(CONCATINATED_FID_FILE, "r"):
        all_fids.add(line.strip("\r\n"))
    all_fids = list(all_fids)
except IOError:
    all_fids = list(set(sub_fids + dlit_fids))
    f = open(CONCATINATED_FID_FILE, "w")
    for fid in all_fids:
        f.write("%s\n" %(fid))
    f.close()

#############
# Filter the subsystem FIDs by organism... we only want OTU genes.
# Unlike the neighborhood analysis, we don't want to include only 
# prokaryotes here.
#############
sys.stderr.write("Filtered list by OTUs...\n")
try:
    otu_fids = readFilteredOtus()
except IOError:
    otu_fids = filterFidsByOtus(all_fids, otus)
    writeFilteredOtus(otu_fids)

#############
# Identify roles for the OTU genes in the organism...
#############
sys.stderr.write("Roles for filtered list...\n")
try:
    otu_fidsToRoles = readFilteredOtuRoles()
except IOError:
    otu_fidsToRoles, otuRolesToFids = fidsToRoles(otu_fids)
    writeFilteredOtuRoles(otu_fidsToRoles)

#############
# Generate a FASTA file
# for the fids in fidsToRoles
#############
sys.stderr.write("Subsystem FASTA file...\n")
try:
    readSubsystemFasta()
except IOError:
    fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys())
    writeSubsystemFasta(fidsToSeqs)

#############
# Get neighborhood info
# for the OTUs (prokaryote only because neighborhoods 
# are generally not conserved for eukaryotes)
#############
#sys.stderr.write("OTU neighborhoods...\n")
#try:
#    fid = open(OTU_NEIGHBORHOOD_FILE, "r")
#    fid.close()
#except IOError:
    # tuplist: [ (contig_id, feature_id, start_location, strand) ]
    # Final file has this and then the roles in a delimited list
#    for prokotu in prokotus:
        # This is mostly because I keep running into incredibly stupid errors.
        # Lets see if I can figure out what the hell is causing them.
#        try:
#            fid = open(os.path.join("data", "OTU", prokotu), "r")
#            fid.close()
#        except IOError:
#            tuplist, fidToRoles = getGenomeNeighborhoodsAndRoles([prokotu])
#            writeOtuNeighborhoods(tuplist, fidToRoles, options.verbose, os.path.join("data", "OTU", prokotu))
#    cmd = "cat %s%s* > %s" %(os.path.join("data", "OTU"), os.sep, OTU_NEIGHBORHOOD_FILE)
#    os.system(cmd)

################
# Complexes --> Roles
# Needed to go from annotation likelihoods
# to reaction likelihoods
# Note that it is easier to go in this direction 
#    Because we need all the roles in a complex to get the probability of that complex.
#
################
sys.stderr.write("Complexes to roles...\n")
try:
    complexToRequiredRoles = readComplexRoles()
except IOError:
    complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(MINN, COUNT)
    writeComplexRoles(complexToRequiredRoles)

########
# reaction --> complex
# Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.
#######

sys.stderr.write("Reactions to complexes...\n")
try:
    rxnToComplexes = readReactionComplex()
except IOError:
    rxnToComplexes, complexesToReactions = reactionComplexLinks(MINN, COUNT)
    writeReactionComplex(rxnToComplexes)

sys.stderr.write("Data gathering done...\n")
