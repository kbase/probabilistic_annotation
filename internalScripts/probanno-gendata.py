#! /usr/bin/python

import sys, time, os, traceback
from biokbase.probabilistic_annotation.DataParser import *
from biokbase.probabilistic_annotation.DataExtractor import *

MINN = 0
COUNT = 50000

def safeRemove(filename):
    try:
        # Check for file existence
        if os.path.exists(filename):
            os.remove(filename)
            
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass
    
def generate_data(config):
    
    # When regenerating the database files, remove all of them first.
    if config["generate_data_option"] == "force":
        sys.stderr.write("Removing all static database files...")
        for filename in DatabaseFiles.values():
            safeRemove(os.path.join(config["data_folder_path"], filename))
        sys.stderr.write("done\n")
    
    sys.stderr.write("Generating static database files in '%s'...\n" %(config["data_folder_path"]))
    
    # Create the data folder if it does not exist.
    if not os.path.exists(config["data_folder_path"]):
        os.makedirs(config["data_folder_path"], 0775)
                
    # Get lists of OTUs
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["otu_id_file"])
    sys.stderr.write("OTU data is in file %s\n  generating file..." %(filename))
    otus, prokotus = getOtuGenomeIds(MINN, COUNT, config)
    writeOtuData(otus, prokotus, config)
    sys.stderr.write("%d otus, %d prokotus...done\n" %(len(otus), len(prokotus)))
    
    # Get a list of subsystem FIDs
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_fid_file"])
    sys.stderr.write("List of subsystem FIDs is in file %s\n  generating file..." %(filename))
    sub_fids = subsystemFids(MINN, COUNT, config)
    writeSubsystemFids(sub_fids, config)
    sys.stderr.write("%d items in list...done\n" %(len(sub_fids)))
    
    # Get a list of Dlit FIDSs
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["dlit_fid_file"])
    sys.stderr.write("List of DLit FIDs is in file %s\n  generating file..." %(filename))
    dlit_fids = getDlitFids(MINN, COUNT, config)
    writeDlitFids(dlit_fids, config)
    sys.stderr.write("%d items in list...done\n" %(len(dlit_fids)))
    
    # Concatenate the two FID lists before filtering
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_file"])
    sys.stderr.write("Combined lists of subsystem and DLit FIDs is in file %s\n  generating file..." %(filename))
    all_fids = list(set(sub_fids + dlit_fids))
    f = open(filename, "w")
    for fid in all_fids:
        f.write("%s\n" %(fid))
    f.close()
    sys.stderr.write("%d items in list...done\n" %(len(all_fids)))
    
    # Identify roles for the OTU genes
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_role_file"])
    sys.stderr.write("List of unfiltered roles is in file %s\n  generating file..." %(filename))
    all_fidsToRoles, all_rolesToFids = fidsToRoles(all_fids, config)
    writeAllFidRoles(all_fidsToRoles, config)
    sys.stderr.write("%d items in list...done\n" %(len(all_fidsToRoles)))
    
    # Filter the subsystem FIDs by organism. We only want OTU genes.
    # Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fid_roles_file"])
    sys.stderr.write("List of filtered OTUs is in file %s\n  generating file..." %(filename))
    otudict = getOtuGenomeDictionary(MINN, COUNT, config)
    otu_fidsToRoles, otu_rolesToFids, missing_roles = filterFidsByOtusBetter(all_fidsToRoles, all_rolesToFids, otudict, config)
    writeFilteredOtuRoles(otu_fidsToRoles, config)
    sys.stderr.write("%d items in list...done\n" %(len(otu_fidsToRoles)))
    
    # Generate a FASTA file for the fids in fidsToRoles
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fasta_file"])
    sys.stderr.write("Subsystem FASTA is in file %s\n  generating file..." %(filename))
    fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys(), config)
    writeSubsystemFasta(fidsToSeqs, config)
    sys.stderr.write("done\n")
    
    # Complexes --> Roles
    # Needed to go from annotation likelihoods
    # to reaction likelihoods
    # Note that it is easier to go in this direction 
    #    Because we need all the roles in a complex to get the probability of that complex.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["complexes_roles_file"])
    sys.stderr.write("Complexes to roles mapping is in file %s\n  generating file..." %(filename))
    complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(MINN, COUNT, config)
    writeComplexRoles(complexToRequiredRoles, config)
    sys.stderr.write("%d items in mapping...done\n" %(len(complexToRequiredRoles)))
    
    # reaction --> complex
    # Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["reaction_complexes_file"])    
    sys.stderr.write("Reactions to complexes mapping is in file %s\n  generating file..." %(filename))
    rxnToComplexes, complexesToReactions = reactionComplexLinks(MINN, COUNT, config)
    writeReactionComplex(rxnToComplexes, config)
    sys.stderr.write("%d items in mapping...done\n" %(len(rxnToComplexes)))
    
    sys.stderr.write("Data gathering done\n")
    return

# Main script function
if __name__ == "__main__":

    # First parameter is the path to the config file.
    if len(sys.argv) != 2:
        if "KB_TOP" not in os.environ:
            sys.stderr.write("KB_TOP environment variable is not set.\n")
            exit(1)
        configFilename = None
    else:
        configFilename = sys.argv[1]
        
    # Read the config from the file.
    config = getConfig(configFilename)
    
    # Generate the static database files.
    try:
        generate_data(config)
        status = "ready"
    except:
        status = "failed"
        sys.stderr.write("\nCaught exception...\n")
        traceback.print_exc(file=sys.stderr)
    
    # Update the status file.
    writeStatusFile(config, status)
    
    exit(0)
