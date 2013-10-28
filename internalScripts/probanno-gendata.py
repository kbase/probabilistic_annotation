#! /usr/bin/python

import sys
import time
import os
import traceback
import argparse
from biokbase.probabilistic_annotation.DataParser import *
from biokbase.probabilistic_annotation.DataExtractor import *

def safeRemove(filename):
    try:
        # Check for file existence
        if os.path.exists(filename):
            os.remove(filename)
            
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass

def timestamp():
    return '%s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

def generate_data(config, force):
    
    # When regenerating the database files, remove all of them first.
    if force:
        sys.stderr.write("Removing all static database files...")
        for filename in DatabaseFiles.values():
            safeRemove(os.path.join(config["data_folder_path"], filename))
        sys.stderr.write("done\n")
    
    sys.stderr.write("Generating static database files in '%s'...\n\n" %(config["data_folder_path"]))
    
    # Create the data folder if it does not exist.
    if not os.path.exists(config["data_folder_path"]):
        os.makedirs(config["data_folder_path"], 0775)
                
    # Get lists of OTUs
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["otu_id_file"])
    sys.stderr.write("OTU data is in file %s\ngenerating file...\n" %(filename))
    otus, prokotus = getOtuGenomeIds(1000, config) # Data build V2 size is 1274
    writeOtuData(otus, prokotus, config)
    sys.stderr.write("%d otus, %d prokotus\ndone at %s\n\n" %(len(otus), len(prokotus), timestamp()))
    
    # Get a list of subsystem FIDs
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_fid_file"])
    sys.stderr.write("List of subsystem FIDs is in file %s\ngenerating file...\n" %(filename))
    sub_fids = subsystemFids(1000, config) # Data build V2 size is 2057
    writeSubsystemFids(sub_fids, config)
    sys.stderr.write("%d items in list\ndone at %s\n\n" %(len(sub_fids), timestamp()))
    
    # Get a list of Dlit FIDSs
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["dlit_fid_file"])
    sys.stderr.write("List of DLit FIDs is in file %s\ngenerating file...\n" %(filename))
    dlit_fids = getDlitFids(5000, config) # Data build V2 size is 12469
    writeDlitFids(dlit_fids, config)
    sys.stderr.write("%d items in list\ndone at %s\n\n" %(len(dlit_fids), timestamp()))
    
    # Concatenate the two FID lists before filtering
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_file"])
    sys.stderr.write("Combined lists of subsystem and DLit FIDs is in file %s\ngenerating file...\n" %(filename))
    all_fids = list(set(sub_fids + dlit_fids))
    f = open(filename, "w")
    for fid in all_fids:
        f.write("%s\n" %(fid))
    f.close()
    sys.stderr.write("%d items in list\ndone at %s\n\n" %(len(all_fids), timestamp()))
    
    # Identify roles for the OTU genes
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_role_file"])
    sys.stderr.write("List of unfiltered roles is in file %s\ngenerating file...\n" %(filename))
    all_fidsToRoles, all_rolesToFids = fidsToRoles(all_fids, config)
    writeAllFidRoles(all_fidsToRoles, config)
    sys.stderr.write("%d items in list\ndone at %s\n\n" %(len(all_fidsToRoles), timestamp()))
    
    # Filter the subsystem FIDs by organism. We only want OTU genes.
    # Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fid_roles_file"])
    sys.stderr.write("List of filtered OTUs is in file %s\ngenerating file...\n" %(filename))
    otudict = getOtuGenomeDictionary(1000, config) # Data build V2 size is 1274
    otu_fidsToRoles, otu_rolesToFids, missing_roles = filterFidsByOtusBetter(all_fidsToRoles, all_rolesToFids, otudict, config)
    writeFilteredOtuRoles(otu_fidsToRoles, config)
    sys.stderr.write("%d items in list\ndone at %s\n\n" %(len(otu_fidsToRoles), timestamp()))
    
    # Generate a FASTA file for the fids in fidsToRoles
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fasta_file"])
    sys.stderr.write("Subsystem FASTA is in file %s\ngenerating file...\n" %(filename))
    fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys(), config)
    writeSubsystemFasta(fidsToSeqs, config)
    sys.stderr.write("done at %s\n\n" %(timestamp()))
    
    # Complexes --> Roles
    # Needed to go from annotation likelihoods
    # to reaction likelihoods
    # Note that it is easier to go in this direction 
    #    Because we need all the roles in a complex to get the probability of that complex.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["complexes_roles_file"])
    sys.stderr.write("Complexes to roles mapping is in file %s\ngenerating file...\n\n" %(filename))
    complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(1000, config) # Data build V2 size is 2369
    writeComplexRoles(complexToRequiredRoles, config)
    sys.stderr.write("%d items in mapping\ndone at %s\n\n" %(len(complexToRequiredRoles), timestamp()))
    
    # reaction --> complex
    # Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["reaction_complexes_file"])    
    sys.stderr.write("Reactions to complexes mapping is in file %s\ngenerating file...\n" %(filename))
    rxnToComplexes, complexesToReactions = reactionComplexLinks(5000, config) # Data build V2 size is 33733
    writeReactionComplex(rxnToComplexes, config)
    sys.stderr.write("%d items in mapping\ndone at %s\n\n" %(len(rxnToComplexes), timestamp()))
    
    sys.stderr.write("Data gathering done\n")
    return

# Main script function
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='probanno-gendata')
    parser.add_argument('configFilename', help='path to config file', action='store', default=None)
    parser.add_argument('--force', help='remove existing static database files first', dest='force', action='store_true', default=False)
    args = parser.parse_args()

    # Read the config from the file.
    config = getConfig(args.configFilename)
    
    # Generate the static database files.
    try:
        generate_data(config, args.force)
        status = "ready"
    except:
        status = "failed"
        sys.stderr.write("\nERROR Caught exception...\n")
        traceback.print_exc(file=sys.stderr)
    
    # Update the status file.
    writeStatusFile(config, status)
    
    exit(0)
