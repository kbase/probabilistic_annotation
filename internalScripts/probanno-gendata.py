#! /usr/bin/python

import sys, time
from biokbase.probabilistic_annotation.DataParser import *
from biokbase.probabilistic_annotation.DataExtractor import *

MINN = 0
COUNT = 50000

def safeRemove(fname, dirname):
    totalfname = os.path.join(dirname, fname)
    try:
        # Check for file existence
        fid = open(totalfname, "r")
        fid.close()
        os.remove(totalfname)
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass
    
def generate_data(config):
    # When regenerating the database files, remove all of them first.
    if config["generate_data_option"] == "force":
        sys.stderr.write("Removing all static database files\n")
        safeRemove(config["otu_id_file"], config["data_folder_path"])
        safeRemove(config["subsystem_fid_file"], config["data_folder_path"])
        safeRemove(config["dlit_fid_file"], config["data_folder_path"])
        safeRemove(config["concatinated_fid_file"], config["data_folder_path"])
        safeRemove(config["subsystem_otu_fid_roles_file"], config["data_folder_path"])
        safeRemove(config["subsystem_otu_fasta_file"], config["data_folder_path"])
        safeRemove(config["subsystem_otu_fasta_file"] + ".psq", config["data_folder_path"]) 
        safeRemove(config["subsystem_otu_fasta_file"] + ".pin", config["data_folder_path"])
        safeRemove(config["subsystem_otu_fasta_file"] + ".phr", config["data_folder_path"])
        safeRemove(config["complexes_roles_file"], config["data_folder_path"])
        safeRemove(config["reaction_complexes_file"], config["data_folder_path"])
    
    sys.stderr.write("Generating static database files...\n")
    
    # Get lists of OTUs
    sys.stderr.write("OTU data...")
    try:
        sys.stderr.write("reading from file...")
        otus, prokotus = readOtuData(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        otus, prokotus = getOtuGenomeIds(MINN, COUNT)
        writeOtuData(otus, prokotus, config)
    sys.stderr.write("done\n")
    
    # Get a list of subsystem FIDs
    sys.stderr.write("List of subsystem FIDS...")
    try:
        sys.stderr.write("reading from file...")
        sub_fids = readSubsystemFids(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        sub_fids = subsystemFids(MINN, COUNT)
        writeSubsystemFids(sub_fids, config)
    sys.stderr.write("done\n")
    
    # Get a list of Dlit FIDSs
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    sys.stderr.write("Getting a list of DLit FIDs...")
    try:
        sys.stderr.write("reading from file...")
        dlit_fids = readDlitFids(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        dlit_fids = getDlitFids(MINN, COUNT)
        writeDlitFids(dlit_fids, config)
    sys.stderr.write("done\n")
    
    # Concatenate the two FID lists before filtering
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    sys.stderr.write("Combining lists of subsystem and DLit FIDS...")
    fn = os.path.join(config["data_folder_path"], config["concatenated_fid_file"])
    try:
        sys.stderr.write("reading from file...")
        all_fids = set()
        for line in open(fn, "r"):
            all_fids.add(line.strip("\r\n"))
        all_fids = list(all_fids)
    except IOError:
        sys.stderr.write("failed...generating file...")
        all_fids = list(set(sub_fids + dlit_fids))
        f = open(fn, "w")
        for fid in all_fids:
            f.write("%s\n" %(fid))
        f.close()
    sys.stderr.write("done\n")
    
    # Identify roles for the OTU genes
    sys.stderr.write("Roles for un-filtered list...")
    try:
        sys.stderr.write("reading from file...")
        all_fidsToRoles, all_rolesToFids = readAllFidRoles(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        all_fidsToRoles, all_rolesToFids = fidsToRoles(all_fids)
        writeAllFidRoles(all_fidsToRoles, config)
    sys.stderr.write("done\n")
    
    # Filter the subsystem FIDs by organism. We only want OTU genes.
    # Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    sys.stderr.write("Filtered list by OTUs...")
    try:
        sys.stderr.write("reading from file...")
        otu_fidsToRoles, otu_rolesToFids = readFilteredOtuRoles(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        otudict = getOtuGenomeDictionary(MINN, COUNT)
        otu_fidsToRoles, otu_rolesToFids, missing_roles = filterFidsByOtusBetter(all_fidsToRoles, all_rolesToFids, otudict)
        writeFilteredOtuRoles(otu_fidsToRoles, config)
    sys.stderr.write("done\n")
    
    # Generate a FASTA file for the fids in fidsToRoles
    sys.stderr.write("Subsystem FASTA file...")
    try:
        sys.stderr.write("reading from file...")
        readSubsystemFasta(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys())
        writeSubsystemFasta(fidsToSeqs, config)
    sys.stderr.write("done\n")
    
    # Complexes --> Roles
    # Needed to go from annotation likelihoods
    # to reaction likelihoods
    # Note that it is easier to go in this direction 
    #    Because we need all the roles in a complex to get the probability of that complex.
    sys.stderr.write("Complexes to roles...")
    try:
        sys.stderr.write("reading from file...")
        complexToRequiredRoles = readComplexRoles(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(MINN, COUNT)
        writeComplexRoles(complexToRequiredRoles, config)
    sys.stderr.write("done\n")
    
    # reaction --> complex
    # Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.    
    sys.stderr.write("Reactions to complexes...")
    try:
        sys.stderr.write("reading from file...")
        rxnToComplexes = readReactionComplex(config)
    except IOError:
        sys.stderr.write("failed...generating file...")
        rxnToComplexes, complexesToReactions = reactionComplexLinks(MINN, COUNT)
        writeReactionComplex(rxnToComplexes, config)
    sys.stderr.write("done\n")
    
    sys.stderr.write("Data gathering done\n")
    return

# Main script function

# First parameter is the path to the config file.
# Read the config from the file.
config = getConfig(sys.argv[1])

# Generate the static database files.
generate_data(config)

# Update the status file.
fid = open(os.path.join(config["data_folder_path"], "status"), "w")
fid.write("ready\ncompleted at %s\n" %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime())))
fid.close()

exit(0)
