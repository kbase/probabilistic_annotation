#! /usr/bin/python

import sys
import time
import os
import traceback
import argparse
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.Helpers import get_config, now, safe_remove
from biokbase.probabilistic_annotation.DataExtractor import CDMExtractor

desc1 = '''
NAME
      pa-gendata-cdm -- generate static database of gene annotations from CDM

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Generate the static database of high-quality gene annotations along with
      files containing intermediate data from KBase central data model.  The
      configFilePath argument specifies the path to the configuration file for
      the service.

      The --force optional argument deletes all existing files before they are
      generated.
'''

desc3 = '''
EXAMPLES
      Generate static database files:
      > pa-gendata gendata.cfg
      
SEE ALSO
      pa-gendata-kegg
      pa-gendata
      pa-savedata
      pa-loaddata

AUTHORS
      Matt Benedict, Mike Mundy 
'''

def generate_data(dataParser, cdmExtractor, force):
    
    print 'Generating CDM static database files in "%s"...' %(config['data_folder_path'])
    print 'Central data model server is at %s' %(config['cdmi_url'])
    print

    # If requested, remove the generated files first.
    if force:
        print 'Removing CDM existing static database files...'
        for filename in dataParser.DataFiles.values():
            safeRemove(filename)
        print 'Done at %s\n' %(now())
    
    # Get list of representative OTU genome IDs.
    print 'Getting list of representative OTU genome IDs at %s' %(now())
    print 'Saving list to file "%s"...' %(dataParser.sources['cdm']['otu_id_file'])
    otus, prokotus = cdmExtractor.getOtuGenomeIds(1000)
    dataParser.writeOtuData(otus, prokotus)
    print 'Found %d OTU genome IDs of which %d are from prokaryotes' %(len(otus), len(prokotus))
    print 'Done at %s\n' %(now())
    del otus, prokotus
    
    # Get a list of subsystem feature IDs (FIDs).
    # Functional annotations from the SEED subsystems are manually curated from
    # multiple sources of information.
    print 'Getting list of subsystem feature IDs at %s' %(now())
    print 'Saving list to file "%s"...' %(dataParser.sources['cdm']['subsystem_fid_file'])
    subsysFids = cdmExtractor.subsystemFids(1000)
    dataParser.writeFeatureIdFile(dataParser.sources['cdm']['subsystem_fid_file'], subsysFids)
    print 'Found %d subsystem feature IDs' %(len(subsysFids))
    print 'Done at %s\n' %(now())
    
    # Get a list of direct literature-supported feature IDs.
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    sys.stderr.write("Getting list of direct literature-supported feature IDs at %s\n" %(now()))
    sys.stderr.write("Saving list to file '%s'\nDownloading from cdmi server...\n" %(dataParser.CdmDataFiles['dlit_fid_file']))
    literatureFids = getDlitFids(5000, config) # Data build V2 size is 12469
    dataParser.writeFeatureIdFile(dataParser.CdmDataFiles['dlit_fid_file'], literatureFids)
    sys.stderr.write("Found %d literature-supported feature IDs\nDone at %s\n\n" %(len(literatureFids), now()))
    
    # Concatenate the two feature ID lists before filtering.
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    sys.stderr.write("Merging lists of subsystem and literature feature IDs at %s\n" %(now()))
    sys.stderr.write("Saving list to file '%s'\nGenerating file...\n" %(dataParser.CdmDataFiles['concatenated_fid_file']))
    allFids = list(set(subsysFids + literatureFids))
    dataParser.writeFeatureIdFile(dataParser.CdmDataFiles['concatenated_fid_file'], allFids)
    sys.stderr.write("Stored %d feature IDs in combined list\nDone at %s\n\n" %(len(allFids), now()))
    del subsysFids, literatureFids
    
    # Identify a role for each feature ID in the concatenated list.
    sys.stderr.write("Getting roles for all feature IDs at %s\n" %(now()))
    sys.stderr.write("Saving mapping of feature ID to roles to file '%s'\nDownloading from cdmi server...\n" %(dataParser.CdmDataFiles['fid_role_file']))
    allFidsToRoles, allRolesToFids = fidsToRoles(allFids, config)
    dataParser.writeFidRoleFile(dataParser.CdmDataFiles['fid_role_file'], allFidsToRoles)
    sys.stderr.write("Stored %d feature ID to roles mappings\nDone at %s\n\n" %(len(allFidsToRoles), now()))
    del allFidsToRoles
    
    # Get a mapping of OTU representative genome IDs to all genomes in the OTU.
    sys.stderr.write("Getting mapping of OTU representative genome IDs to all genomes in OTU at %s\n" %(now()))
    sys.stderr.write("Downloading from cdmi server...\n")
    otuGenomes = getOtuGenomeDictionary(1000, config) # Data build V2 size is 1274
    sys.stderr.write("Found %d representative OTU genome IDs\nDone at %s\n\n" %(len(otuGenomes), now()))
    
    # Filter the feature IDs by organism. We only want one feature ID from each OTU for each
    # functional role.  Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    sys.stderr.write("Filtering list of feature IDs so there is one protein from each OTU for each functional role at %s\n" %(now()))
    sys.stderr.write("Saving list of filtered feature IDs in file '%s'\nQuerying cdmi server...\n" %(dataParser.CdmDataFiles['otu_fid_role_file']))
    otuFidsToRoles, otuRolesToFids = filterFidsByOtusOptimized(allFids, allRolesToFids, otuGenomes, config)
    dataParser.writeFidRoleFile(dataParser.CdmDataFiles['otu_fid_role_file'], otuFidsToRoles)
    sys.stderr.write("Stored %d feature ID to role mappings\nDone at %s\n\n" %(len(otuFidsToRoles), now()))
    del allFids, otuRolesToFids, otuGenomes
    
    # Generate a FASTA file for the feature IDs in filtered list and make a BLAST database.
    sys.stderr.write("Getting amino acid sequences for filtered feature IDs at %s\n" %(now()))
    sys.stderr.write("Downloading from cdmi server...\n")
    fidsToSeqs = fidsToSequences(otuFidsToRoles.keys(), config)
    sys.stderr.write("Writing amino acid sequences to FASTA file '%s'\nGenerating file and making search database...\n" %(dataParser.DataFiles['protein_fasta_file']))
    dataParser.writeProteinFastaFile(fidsToSeqs)
#    dataParser.buildSearchDatabase()
    sys.stderr.write("Done at %s\n\n" %(now()))
    del otuFidsToRoles, fidsToSeqs
    
    # Create a mapping of complexes to roles which is needed to go from annotation likelihoods to
    # reaction likelihoods.  Note that it is easier to go in this direction because we need all
    # the roles in a complex to get the probability of that complex.
    sys.stderr.write("Getting mapping of complex to roles at %s\n" %(now()))
    sys.stderr.write("Saving complex to roles mapping in file '%s'\nDownloading from cdmi server...\n" %(dataParser.CdmDataFiles['complex_role_file']))
    complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(1000, config) # Data build V2 size is 2369
    dataParser.writeComplexRoleFile(dataParser.CdmDataFiles['complex_role_file'], complexToRequiredRoles)
    sys.stderr.write("Stored %d complex to roles mappings\nDone at %s\n\n" %(len(complexToRequiredRoles), now()))
    del complexToRequiredRoles, requiredRolesToComplexes
    
    # Create a mapping of reactions to complexes.  Note that it is easier to go in this direction since
    # we'll be filtering multiple complexes down to a single reaction.
    sys.stderr.write("Getting mapping of reaction to complexes at %s\n" %(now()))
    sys.stderr.write("Saving reaction to complexes mapping in file '%s'\nDownloading from cdmi server...\n" %(dataParser.CdmDataFiles['reaction_complex_file']))
    reactionToComplexes, complexesToReactions = reactionComplexLinks(5000, config) # Data build V2 size is 33733
    dataParser.writeReactionComplexFile(dataParser.CdmDataFiles['reaction_complex_file'], reactionToComplexes)
    sys.stderr.write("Stored %d reaction to complexes mappings\nDone at %s\n\n" %(len(reactionToComplexes), now()))
    del reactionToComplexes, complexesToReactions
    
    sys.stderr.write("Done generating static database files\n")
    return

# Main script function
if __name__ == '__main__':
    # Force output to be unbuffered.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)

    # Parse arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-gendata', epilog=desc3)
    parser.add_argument('configFilePath', help='path to configuration file', action='store', default=None)
    parser.add_argument('--force', help='remove existing static database files first', dest='force', action='store_true', default=False)
    parser.add_argument('--makedb', help='only make the protein search database', dest='makeDB', action='store_true', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)

    # Create a DataParser object for working with the static database files (the
    # data folder is created if it does not exist).
    dataParser = DataParser(config)

    # Create a CDMExtractor object for getting data from the central data model.
    cdmExtractor = CDMExtractor(config)

    # Update the status file.
    dataParser.writeStatusFile('building')

    # Generate the static database files.
    try:
        generate_data(dataParser, cdmExtractor, args.force)
        status = "ready"
    except:
        status = "failed"
        sys.stderr.write("\nERROR Caught exception...\n")
        traceback.print_exc(file=sys.stderr)
    
    # Update the status file.
    dataParser.writeStatusFile(status)
    
    exit(0)
