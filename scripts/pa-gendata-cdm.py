#! /usr/bin/python

import sys
import time
import os
import traceback
import argparse
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.Helpers import get_config, now, safe_remove
from biokbase.probabilistic_annotation.CDMExtractor import CDMExtractor

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

      The --step optional argument starts the process at the specified step,
      which can be used to restart if there is an error.

      The method to generate intermediate data from the central data model is as
      follows: (1) Get a list of representative OTU genome IDs. (2) Get a list
      of subsystem feature IDs. (3) Get a list of direct literature-supported
      feature IDs. (4) Merge the two lists of feature IDs. (5) Create a mapping
      from feature IDs to roles. (6) Filter the feature IDs so there is one
      protein from each OTU for each functional role. (7) Get the amino acid
      sequence for each filtered feature ID. (8) Create a mapping of complexes
      to roles. (9) Create a mapping of reactions to complexes.
'''

desc3 = '''
EXAMPLES
      Generate static database files:
      > pa-gendata-cdm gendata.cfg
      
SEE ALSO
      pa-gendata-kegg
      pa-gendata
      pa-savedata
      pa-loaddata

AUTHORS
      Matt Benedict, Mike Mundy 
'''

def generate_data(dataParser, cdmExtractor, force, startStep):
    ''' Generate intermediate database files from central data model.

        @param dataParser: DataParser object for reading and writing data files
        @param cdmExtractor: CDMExtractor object for extracting data from CDM
        @param force: When true, remove existing files first
        @param startStep: Start generating data from this step
        @return Nothing
    '''
    
    # If requested, remove the generated files first.
    if force:
        print 'Removing existing CDM static database files...'
        for filename in dataParser.sources['cdm'].values():
            safe_remove(filename)
        print 'Done at %s\n' %(now())
    
    # Get list of representative OTU genome IDs.
    step = 1
    print 'Getting list of representative OTU genome IDs at %s' %(now())
    if step >= startStep:
        print 'Saving list to file "%s"...' %(dataParser.sources['cdm']['otu_id_file'])
        otus, prokotus = cdmExtractor.getOtuGenomeIds(1000)
        dataParser.writeOtuData(dataParser.sources['cdm']['otu_id_file'], otus, prokotus)
        print 'Stored %d OTU genome IDs of which %d are from prokaryotes' %(len(otus), len(prokotus))
    else:
        print 'Loading list from file "%s"...' %(dataParser.sources['cdm']['otu_id_file'])
        otus, prokotus = dataParser.readOtuData(dataParser.sources['cdm']['otu_id_file'])
        print 'Retrieved %d OTU genome IDs of which %d are from prokaryotes' %(len(otus), len(prokotus))
    del otus, prokotus
    print 'Done at %s\n' %(now())
    step += 1

    # Get a list of subsystem feature IDs (FIDs).
    # Functional annotations from the SEED subsystems are manually curated from
    # multiple sources of information.
    print 'Getting list of subsystem feature IDs at %s' %(now())
    if step >= startStep:
        print 'Saving list to file "%s"...' %(dataParser.sources['cdm']['subsystem_fid_file'])
        subsysFids = cdmExtractor.getSubsystemFids(1000)
        dataParser.writeFeatureIdFile(dataParser.sources['cdm']['subsystem_fid_file'], subsysFids)
        print 'Stored %d subsystem feature IDs' %(len(subsysFids))
    else:
        print 'Loading list from file "%s"...' %(dataParser.sources['cdm']['subsystem_fid_file'])
        subsysFids = dataParser.readFeatureIdFile(dataParser.sources['cdm']['subsystem_fid_file'])
        print 'Retrieved %d subsystem feature IDs' %(len(subsysFids))
    print 'Done at %s\n' %(now())
    step += 1

    # Get a list of direct literature-supported feature IDs.
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    print 'Getting list of direct literature-supported feature IDs at %s' %(now())
    if step >= startStep:
        print 'Saving list to file "%s"...' %(dataParser.sources['cdm']['dlit_fid_file'])
        literatureFids = cdmExtractor.getDlitFids(1000)
        dataParser.writeFeatureIdFile(dataParser.sources['cdm']['dlit_fid_file'], literatureFids)
        print 'Stored %d literature-supported feature IDs' %(len(literatureFids))
    else:
        print 'Loading list from file "%s"...' %(dataParser.sources['cdm']['dlit_fid_file'])
        literatureFids = dataParser.readFeatureIdFile(dataParser.sources['cdm']['dlit_fid_file'])
        print 'Retrieved %d literature-supported feature IDs' %(len(literatureFids))
    print 'Done at %s\n' %(now())
    step += 1

    # Concatenate the two feature ID lists before filtering.
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    print 'Merging lists of subsystem and literature feature IDs at %s' %(now())
    if step >= startStep:
        print 'Saving list to file "%s"...' %(dataParser.sources['cdm']['concatenated_fid_file'])
        allFids = list(set(subsysFids + literatureFids))
        dataParser.writeFeatureIdFile(dataParser.sources['cdm']['concatenated_fid_file'], allFids)
        print 'Stored %d feature IDs in merged list' %(len(allFids))
    else:
        print 'Loading list from file "%s"...' %(dataParser.sources['cdm']['concatenated_fid_file'])
        allFids = dataParser.readFeatureIdFile(dataParser.sources['cdm']['concatenated_fid_file'])
        print 'Retrieved %d feature IDs in merged list' %(len(allFids))
    print 'Done at %s\n' %(now())
    del subsysFids, literatureFids
    step += 1

    # Identify a role for each feature ID in the concatenated list.
    print 'Getting roles for all feature IDs at %s' %(now())
    if step >= startStep:
        print 'Saving mapping of feature ID to roles to file "%s"...' %(dataParser.sources['cdm']['fid_role_file'])
        allFidsToRoles, allRolesToFids = cdmExtractor.mapFidsToRoles(allFids)
        dataParser.writeFidRoleFile(dataParser.sources['cdm']['fid_role_file'], allFidsToRoles)
        print 'Stored %d mappings of feature ID to roles' %(len(allFidsToRoles))
    else:
        print 'Loading mapping of feature ID to roles from file "%s"...' %(dataParser.sources['cdm']['fid_role_file'])
        allFidsToRoles, allRolesToFids = dataParser.readFidRoleFile(dataParser.sources['cdm']['fid_role_file'])
        print 'Retrieved %d mappings of feature ID to roles' %(len(allFidsToRoles))
    print 'Done at %s\n' %(now())
    del allFidsToRoles
    step += 1

    # Get a mapping of OTU representative genome IDs to all genomes in the OTU.
    print 'Getting mapping of OTU representative genome IDs to all genomes in OTU at %s' %(now())
    otuGenomes = cdmExtractor.getOtuGenomeDictionary(1000)
    print 'Found %d representative OTU genome IDs' %(len(otuGenomes))
    print 'Done at %s\n' %(now())

    # Filter the feature IDs by organism. We only want one feature ID from each OTU for each
    # functional role.  Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    print 'Filtering list of feature IDs so there is one protein from each OTU for each functional role at %s' %(now())
    if step >= startStep:
        print 'Saving list of filtered feature IDs in file "%s"...' %(dataParser.sources['cdm']['otu_fid_role_file'])
        otuFidsToRoles, otuRolesToFids = cdmExtractor.filterFidsByOtusOptimized(allFids, allRolesToFids, otuGenomes)
        dataParser.writeFidRoleFile(dataParser.sources['cdm']['otu_fid_role_file'], otuFidsToRoles)
        print 'Stored %d feature ID to role mappings' %(len(otuFidsToRoles))
    else:
        print 'Loading list of filtered feature IDs from file "%s"...' %(dataParser.sources['cdm']['otu_fid_role_file'])
        otuFidsToRoles, otuRolesToFids = dataParser.readFidRoleFile(dataParser.sources['cdm']['otu_fid_role_file'])
        print 'Retrieved %d feature ID to role mappings' %(len(otuFidsToRoles))
    print 'Done at %s\n' %(now())
    del allFids, otuRolesToFids, otuGenomes
    step += 1

    # Generate a FASTA file for the feature IDs in filtered list.
    print 'Getting amino acid sequences for filtered feature IDs at %s' %(now())
    if step >= startStep:
        print 'Saving amino acid sequences to FASTA file "%s"...' %(dataParser.sources['cdm']['protein_fasta_file'])
        fidsToSeqs = cdmExtractor.getAminoAcidSequences(otuFidsToRoles.keys())
        dataParser.writeProteinFastaFile(dataParser.sources['cdm']['protein_fasta_file'], fidsToSeqs)
        print 'Stored %d amino acid sequences' %(len(fidsToSeqs))
        del fidsToSeqs
    print 'Done at %s\n' %(now())
    del otuFidsToRoles
    step += 1

    # Create a mapping of complexes to roles which is needed to go from annotation likelihoods to
    # reaction likelihoods.  Note that it is easier to go in this direction because we need all
    # the roles in a complex to get the probability of that complex.
    print 'Getting mapping of complex to role at %s' %(now())
    if step >= startStep:
        print 'Saving complex to role mapping in file "%s"...' %(dataParser.sources['cdm']['complex_role_file'])
        complexToRequiredRoles, requiredRolesToComplexes = cdmExtractor.mapComplexToRole(1000)
        dataParser.writeComplexRoleFile(dataParser.sources['cdm']['complex_role_file'], complexToRequiredRoles)
        print 'Stored %d complex to role mappings' %(len(complexToRequiredRoles))
        del complexToRequiredRoles, requiredRolesToComplexes
    print 'Done at %s\n' %(now())
    step += 1

    # Create a mapping of reactions to complexes.  Note that it is easier to go in this direction since
    # we'll be filtering multiple complexes down to a single reaction.
    print 'Getting mapping of reaction to complex at %s' %(now())
    if step >= startStep:
        print 'Saving reaction to complex mapping in file "%s"...' %(dataParser.sources['cdm']['reaction_complex_file'])
        reactionToComplexes, complexesToReactions = cdmExtractor.mapReactionToComplex(5000)
        dataParser.writeReactionComplexFile(dataParser.sources['cdm']['reaction_complex_file'], reactionToComplexes)
        print 'Stored %d reaction to complex mappings' %(len(reactionToComplexes))
        del reactionToComplexes, complexesToReactions
    print 'Done at %s\n' %(now())
    step += 1

    print 'Done generating CDM static database files'
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
    parser.add_argument('--step', help='start generating database files from this step', dest='step', action='store', type=int, default=1)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)
    print 'Generating CDM static database files in directory "%s"' %(config['data_folder_path'])
    print 'Central data model server is at %s' %(config['cdmi_url'])
    print

    # Create a DataParser object for working with the static database files (the
    # data folder is created if it does not exist).
    dataParser = DataParser(config)

    # Create a CDMExtractor object for getting data from the central data model.
    cdmExtractor = CDMExtractor(config)

    # Generate the intermediate static database files from the central data model.
    # Note the status file is not used when generating intermediate files.
    try:
        generate_data(dataParser, cdmExtractor, args.force, args.step)
    except:
        sys.stderr.write('\nERROR Caught exception at %s...\n' %(now()))
        traceback.print_exc(file=sys.stderr)
    
    exit(0)
