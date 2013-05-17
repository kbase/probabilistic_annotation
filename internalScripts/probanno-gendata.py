#! /usr/bin/python

import sys, time, os, traceback
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
    sys.stderr.write("OTU data is in file %s\n" %(os.path.join(config["data_folder_path"], config["otu_id_file"])))
    try:
        sys.stderr.write("  reading from file...")
        otus, prokotus = readOtuData(config)
        sys.stderr.write("%d otus, %d prokotus..." %(len(otus), len(prokotus)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        otus, prokotus = getOtuGenomeIds(MINN, COUNT, config)
        writeOtuData(otus, prokotus, config)
        sys.stderr.write("%d otus, %d prokotus..." %(len(otus), len(prokotus)))
    sys.stderr.write("done\n")
    
    # Get a list of subsystem FIDs
    sys.stderr.write("List of subsystem FIDs is in file %s\n" %(os.path.join(config["data_folder_path"], config["subsystem_fid_file"])))
    try:
        sys.stderr.write("  reading from file...")
        sub_fids = readSubsystemFids(config)
        sys.stderr.write("%d items in list..." %(len(sub_fids)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        sub_fids = subsystemFids(MINN, COUNT, config)
        writeSubsystemFids(sub_fids, config)
        sys.stderr.write("%d items in list..." %(len(sub_fids)))
    sys.stderr.write("done\n")
    
    # Get a list of Dlit FIDSs
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    sys.stderr.write("List of DLit FIDs is in file %s\n" %(os.path.join(config["data_folder_path"], config["dlit_fid_file"])))
    try:
        sys.stderr.write("  reading from file...")
        dlit_fids = readDlitFids(config)
        sys.stderr.write("%d items in list..." %(len(dlit_fids)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        dlit_fids = getDlitFids(MINN, COUNT, config)
        writeDlitFids(dlit_fids, config)
        sys.stderr.write("%d items in list..." %(len(dlit_fids)))
    sys.stderr.write("done\n")
    
    # Concatenate the two FID lists before filtering
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    fn = os.path.join(config["data_folder_path"], config["concatenated_fid_file"])
    sys.stderr.write("Combined lists of subsystem and DLit FIDs is in file %s\n" %(fn))
    try:
        sys.stderr.write("  reading from file...")
        all_fids = set()
        for line in open(fn, "r"):
            all_fids.add(line.strip("\r\n"))
        all_fids = list(all_fids)
        sys.stderr.write("%d items in list..." %(len(all_fids)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        all_fids = list(set(sub_fids + dlit_fids))
        f = open(fn, "w")
        for fid in all_fids:
            f.write("%s\n" %(fid))
        f.close()
    sys.stderr.write("done\n")
    
    # Identify roles for the OTU genes
    sys.stderr.write("List of unfiltered roles is in file %s\n" %(os.path.join(config["data_folder_path"], config["concatenated_fid_role_file"])))
    try:
        sys.stderr.write("  reading from file...")
        all_fidsToRoles, all_rolesToFids = readAllFidRoles(config)
        sys.stderr.write("%d items in list..." %(len(all_fidsToRoles)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        all_fidsToRoles, all_rolesToFids = fidsToRoles(all_fids, config)
        writeAllFidRoles(all_fidsToRoles, config)
        sys.stderr.write("%d items in list..." %(len(all_fidsToRoles)))
    sys.stderr.write("done\n")
    
    # Filter the subsystem FIDs by organism. We only want OTU genes.
    # Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    sys.stderr.write("List of filtered OTUs is in file %s\n" %(os.path.join(config["data_folder_path"], config["subsystem_otu_fid_roles_file"])))
    try:
        sys.stderr.write("  reading from file...")
        otu_fidsToRoles, otu_rolesToFids = readFilteredOtuRoles(config)
        sys.stderr.write("%d items in list..." %(len(otu_fidsToRoles)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        otudict = getOtuGenomeDictionary(MINN, COUNT, config)
        otu_fidsToRoles, otu_rolesToFids, missing_roles = filterFidsByOtusBetter(all_fidsToRoles, all_rolesToFids, otudict, config)
        writeFilteredOtuRoles(otu_fidsToRoles, config)
        sys.stderr.write("%d items in list..." %(len(otu_fidsToRoles)))
    sys.stderr.write("done\n")
    
    # Generate a FASTA file for the fids in fidsToRoles
    sys.stderr.write("Subsystem FASTA is in file %s\n" %(os.path.join(config["data_folder_path"], config["subsystem_otu_fasta_file"])))
    try:
        sys.stderr.write("  reading from file...")
        readSubsystemFasta(config)
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys(), config)
        writeSubsystemFasta(fidsToSeqs, config)
    sys.stderr.write("done\n")
    
    # Complexes --> Roles
    # Needed to go from annotation likelihoods
    # to reaction likelihoods
    # Note that it is easier to go in this direction 
    #    Because we need all the roles in a complex to get the probability of that complex.
    sys.stderr.write("Complexes to roles mapping is in file %s\n" %(os.path.join(config["data_folder_path"], config["complexes_roles_file"])))
    try:
        sys.stderr.write("  reading from file...")
        complexToRequiredRoles = readComplexRoles(config)
        sys.stderr.write("%d items in mapping..." %(len(complexToRequiredRoles)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(MINN, COUNT, config)
        writeComplexRoles(complexToRequiredRoles, config)
        sys.stderr.write("%d items in mapping..." %(len(complexToRequiredRoles)))
    sys.stderr.write("done\n")
    
    # reaction --> complex
    # Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.    
    sys.stderr.write("Reactions to complexes mapping is in file %s\n" %(os.path.join(config["data_folder_path"], config["reaction_complexes_file"])))
    try:
        sys.stderr.write("  reading from file...")
        rxnToComplexes = readReactionComplex(config)
        sys.stderr.write("%d items in mapping..." %(len(rxnToComplexes)))
    except IOError:
        sys.stderr.write("failed\n  generating file...")
        rxnToComplexes, complexesToReactions = reactionComplexLinks(MINN, COUNT, config)
        writeReactionComplex(rxnToComplexes, config)
        sys.stderr.write("%d items in mapping..." %(len(rxnToComplexes)))
    sys.stderr.write("done\n")
    
    sys.stderr.write("Data gathering done\n")
    return

# Main script function

# First parameter is the path to the config file.
# Read the config from the file.
config = getConfig(sys.argv[1])

# Generate the static database files.
try:
    generate_data(config)
    status = "ready"
except:
    status = "failed"
    sys.stderr.write("\nCaught exception...\n")
    traceback.print_exc(file=sys.stderr)

# Update the status file.
fid = open(os.path.join(config["data_folder_path"], config["status_file"]), "w")
fid.write("%s\ncompleted at %s\n" %(status, time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime())))
fid.close()

exit(0)
