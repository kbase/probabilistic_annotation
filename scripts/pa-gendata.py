#! /usr/bin/python

import sys
import time
import os
import traceback
import argparse
from biokbase.probabilistic_annotation.DataParser import *
from biokbase.probabilistic_annotation.DataExtractor import *

desc1 = '''
NAME
      pa-gendata -- generate static database of gene annotations

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Generate the static database of high-quality gene annotations along with
      files containing intermediate data.  The dataFolderPath argument specifies
      the path to the directory where the database and related files are stored.

      The --force optional argument deletes all existing files before they are
      generated.
      
      The --separator optional argument specifies the string to use as a
      separator for list items.  The separator string must not occur in any of
      the values of the list items.

      The --cdmi-url optional argument specifies an alternate URL for the
      central data model interface service.
'''

desc3 = '''
EXAMPLES
      Generate static database files:
      > pa-gendata /kb/deployment/services/probabilistic_annotation/data
      
SEE ALSO
      pa-savedata

AUTHORS
      Matt Benedict, Mike Mundy 
'''

def safeRemove(filename):
    try:
        # Check for file existence
        if os.path.exists(filename):
            os.remove(filename)
            
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass

def generate_data(config, force):
    
    # Create the data folder if it does not exist.
    if not os.path.exists(config["data_folder_path"]):
        os.makedirs(config["data_folder_path"], 0775)
                
    # Update the status file.
    writeStatusFile(config, 'building')

    # When regenerating the database files, remove all of them first.
    if force:
        sys.stderr.write("Removing all static database files...")
        for filename in DatabaseFiles.values():
            safeRemove(os.path.join(config["data_folder_path"], filename))
        sys.stderr.write("done\n")
    
    sys.stderr.write("Generating static database files in '%s'...\n" %(config["data_folder_path"]))
    sys.stderr.write("Central data model server is at %s\n\n" %(config['cdmi_url']))
    
    # Get list of representative OTU genome IDs.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["otu_id_file"])
    sys.stderr.write("Getting list of representative OTU genome IDs at %s\n" %(timestamp()))
    sys.stderr.write("Saving list to file '%s'\nDownloading from cdmi server...\n" %(filename))
    otus, prokotus = getOtuGenomeIds(1000, config) # Data build V2 size is 1274
    writeOtuData(otus, prokotus, config)
    sys.stderr.write("Found %d OTU genome IDs of which %d are from prokaryotes\nDone at %s\n\n" %(len(otus), len(prokotus), timestamp()))
    del otus, prokotus
    
    # Get a list of subsystem feature IDs (FIDs).
    # Functional annotations from the SEED subsystems are manually curated from
    # multiple sources of information.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_fid_file"])
    sys.stderr.write("Getting list of subsystem feature IDs at %s\n" %(timestamp()))
    sys.stderr.write("Saving list to file '%s'\nDownloading from cdmi server...\n" %(filename))
    subsysFids = subsystemFids(1000, config) # Data build V2 size is 2057
    writeSubsystemFids(subsysFids, config)
    sys.stderr.write("Found %d subsystem feature IDs\nDone at %s\n\n" %(len(subsysFids), timestamp()))
    
    # Get a list of direct literature-supported feature IDs.
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["dlit_fid_file"])
    sys.stderr.write("Getting list of direct literature-supported feature IDs at %s\n" %(timestamp()))
    sys.stderr.write("Saving list to file '%s'\nDownloading from cdmi server...\n" %(filename))
    literatureFids = getDlitFids(5000, config) # Data build V2 size is 12469
    writeDlitFids(literatureFids, config)
    sys.stderr.write("Found %d literature-supported feature IDs\nDone at %s\n\n" %(len(literatureFids), timestamp()))
    
    # Concatenate the two feature ID lists before filtering.
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_file"])
    sys.stderr.write("Merging lists of subsystem and literature feature IDs at %s\n" %(timestamp()))
    sys.stderr.write("Saving list to file '%s'\nGenerating file...\n" %(filename))
    allFids = list(set(subsysFids + literatureFids))
    f = open(filename, "w")
    for fid in allFids:
        f.write("%s\n" %(fid))
    f.close()
    sys.stderr.write("Stored %d feature IDs in combined list\nDone at %s\n\n" %(len(allFids), timestamp()))
    del subsysFids, literatureFids
    
    # Identify a role for each feature ID in the concatenated list.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_role_file"])
    sys.stderr.write("Getting roles for all feature IDs at %s\n" %(timestamp()))
    sys.stderr.write("Saving mapping of feature ID to roles to file '%s'\nDownloading from cdmi server...\n" %(filename))
    allFidsToRoles, allRolesToFids = fidsToRoles(allFids, config)
    writeAllFidRoles(allFidsToRoles, config)
    sys.stderr.write("Stored %d feature ID to roles mappings\nDone at %s\n\n" %(len(allFidsToRoles), timestamp()))
    del allFidsToRoles
    
    # Get a mapping of OTU representative genome IDs to all genomes in the OTU.
    sys.stderr.write("Getting mapping of OTU representative genome IDs to all genomes in OTU at %s\n" %(timestamp()))
    sys.stderr.write("Downloading from cdmi server...\n")
    otuGenomes = getOtuGenomeDictionary(1000, config) # Data build V2 size is 1274
    sys.stderr.write("Found %d representative OTU genome IDs\nDone at %s\n\n" %(len(otuGenomes), timestamp()))
    
    # Filter the feature IDs by organism. We only want one feature ID from each OTU for each
    # functional role.  Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fid_roles_file"])
    sys.stderr.write("Filtering list of feature IDs so there is one protein from each OTU for each functional role at %s\n" %(timestamp()))
    sys.stderr.write("Saving list of filtered feature IDs in file '%s'\nQuerying cdmi server...\n" %(filename))
    otuFidsToRoles, otuRolesToFids = filterFidsByOtusOptimized(allFids, allRolesToFids, otuGenomes, config)
    writeFilteredOtuRoles(otuFidsToRoles, config)
    sys.stderr.write("Stored %d feature ID to role mappings\nDone at %s\n\n" %(len(otuFidsToRoles), timestamp()))
    del allFids, otuRolesToFids, otuGenomes
    
    # Generate a FASTA file for the feature IDs in filtered list and make a BLAST database.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fasta_file"])
    sys.stderr.write("Getting amino acid sequences for filtered feature IDs at %s\n" %(timestamp()))
    sys.stderr.write("Downloading from cdmi server...\n")
    fidsToSeqs = fidsToSequences(otuFidsToRoles.keys(), config)
    sys.stderr.write("Writing amino acid sequences to FASTA file '%s'\nGenerating file and making BLAST database...\n" %(filename))
    writeSubsystemFasta(fidsToSeqs, config)
    sys.stderr.write("Done at %s\n\n" %(timestamp()))
    del otuFidsToRoles, fidsToSeqs
    
    # Create a mapping of complexes to roles which is needed to go from annotation likelihoods to
    # reaction likelihoods.  Note that it is easier to go in this direction because we need all
    # the roles in a complex to get the probability of that complex.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["complexes_roles_file"])
    sys.stderr.write("Getting mapping of complex to roles at %s\n" %(timestamp()))
    sys.stderr.write("Saving complex to roles mapping in file '%s'\nDownloading from cdmi server...\n" %(filename))
    complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(1000, config) # Data build V2 size is 2369
    writeComplexRoles(complexToRequiredRoles, config)
    sys.stderr.write("Stored %d complex to roles mappings\nDone at %s\n\n" %(len(complexToRequiredRoles), timestamp()))
    del complexToRequiredRoles, requiredRolesToComplexes
    
    # Create a mapping of reactions to complexes.  Note that it is easier to go in this direction since
    # we'll be filtering multiple complexes down to a single reaction.
    filename = os.path.join(config["data_folder_path"], DatabaseFiles["reaction_complexes_file"])
    sys.stderr.write("Getting mapping of reaction to complexes at %s\n" %(timestamp()))    
    sys.stderr.write("Saving reaction to complexes mapping in file '%s'\nDownloading from cdmi server...\n" %(filename))
    reactionToComplexes, complexesToReactions = reactionComplexLinks(5000, config) # Data build V2 size is 33733
    writeReactionComplex(reactionToComplexes, config)
    sys.stderr.write("Stored %d reaction to complexes mappings\nDone at %s\n\n" %(len(reactionToComplexes), timestamp()))
    del reactionToComplexes, complexesToReactions
    
    sys.stderr.write("Done generating static database files\n")
    return

# Main script function
if __name__ == "__main__":

    # Parse arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-gendata', epilog=desc3)
    parser.add_argument('dataFolderPath', help='path to directory for storing static database files', action='store', default=None)
    parser.add_argument('--cdmi-url', help='url for central data model interface service', action='store', dest='cdmiURL', default='http://kbase.us/services/cdmi_api/')
    parser.add_argument('--force', help='remove existing static database files first', dest='force', action='store_true', default=False)
    parser.add_argument('--separator', help='character string not found in any roles and used to split lists', dest='separator', action='store', default='///')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Build a configuration dictionary.
    config = dict()
    config['data_folder_path'] = args.dataFolderPath
    config['cdmi_url'] = args.cdmiURL
    config['separator'] = args.separator

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
