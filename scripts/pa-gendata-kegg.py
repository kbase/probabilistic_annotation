#! /usr/bin/python

import argparse
import os
import sys
import traceback
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.kegg.KEGGExtractor import KEGGExtractor
from biokbase.probabilistic_annotation.Helpers import get_config, now

desc1 = '''
NAME
      pa-gendata-kegg -- generate static database of gene annotations from KEGG

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Generate the static database of high-quality gene annotations along with
      files containing intermediate data from KEGG.  The configFilePath argument
      specifies the path to the configuration file for the service.  The
      rxnfilename argument specifies the path to a file with a list of reactions.
      The reactions must be defined in a Biochemistry object for the RxnProbs
      objects created by the probabilistic annotation service to be usable by
      the fba modeling service.

      The --source optional argument specifies the type of the source reaction
      file.  The supported types are "template" for a reaction file created by
      the fba-export-template command and "modelseed" for a master reaction file
      from the ModelSEEDDatabase repository.

      The --biochem optional argument specifies the ID of a Biochemistry object
      for looking up KEGG reaction IDs.  The --biochemws optional argument specifies
      the workspace of the Biochemistry object.  The --increment optional argument
      specifies the size of the increment when retrieving lists from the fba
      modeling service.

      The --reactiondb optional argument specifies the path to a KEGG reaction
      database.  The default is "reaction.db" in the current working directory.
      The --enzymedb optional argument specifies the path to a KEGG enzyme
      database.  The default is "enzyme.db" in the current working directory.
      The --organismdb optional argument specifies the path to a KEGG organism
      database.  The default is "organism.db" in the current working directory.
      When the --no-update optional argument is specified, missing reactions or
      enzymes are not retrieved from the KEGG web service.

      The --force optional argument deletes all existing files before they are
      generated.

      The --step optional argument starts the process at the specified step,
      which can be used to restart if there is an error.

      The method to generate intermediate data from KEGG is as follows: (1) Lookup
      reactions in input reaction file in the Biochemistry object and find all
      of the reactions that have a KEGG reaction ID. (2) If requested, update
      the local KEGG reaction and enzyme databases with missing values. (3) Create
      a mapping from reactions to complexes.  Each enzyme listed in a reaction
      is assigned a KBase complex ID starting from kb|cpx.10000. (4) Create a
      mapping from complexes to roles.  Each enzyme assigned a complex ID is
      assigned a role using the enzyme name and EC number. (5) For each enzyme
      assigned a role, get a list of features that are from prokaryote organisms
      that are OTU representatives.  Use the KEGG web service to get the amino
      acid sequence for the feature.

      Note that you must have local copies of the KEGG databases.
'''

desc3 = '''
EXAMPLES
      > pa-gendata-kegg gendata.cfg template.tsv

SEE ALSO
      pa-gendata-cdm
      pa-gendata
      pa-savedata
      pa-loaddata

AUTHORS
      Mike Mundy 
'''

def generate_data(dataParser, keggExtractor, args):
    ''' Generate intermediate database files from KEGG.

        @param dataParser: DataParser object for reading and writing data files
        @param keggExtractor: KEGGExtractor object for extracting data from KEGG
        @param args: Input arguments to script
        @return Nothing
    '''

    # If requested, remove the generated files first.
    if args.force:
        print 'Removing existing KEGG static database files...'
        for filename in dataParser.sources['kegg'].values():
            safe_remove(filename)
        print 'Done at %s\n' %(now())

    # First step is map KBase reaction IDs to KEGG reaction IDs.
    step = 1
    print 'Getting mapping of KBase reaction IDs to KEGG reaction IDs at %s' %(now())
    if step >= args.step:
        print 'Saving mapping to file "%s"...' %(dataParser.sources['kegg']['reaction_file'])
        kbaseRxnIdList = list()
    
        if args.source == 'template':
            # Read the input Model Template file.  Build a list of reactions from the cytosol
            # and extracellular compartments. In the Model Template text file, the first field
            # contains the KBase reaction ID.
            with open(args.rxnfilename, 'r') as handle:
                header = handle.readline() # Throw away the header line
                for line in handle:
                    fields = line.strip('\n').split('\t')
                    if (fields[1] == 'c' or fields[1] == 'e'):
                        kbaseRxnIdList.append(fields[0])
            print 'Found %d KBase reactions in ModelTemplate' %(len(kbaseRxnIdList))
    
        elif args.source == 'modelseed':
            # Read the input master reactions file from ModelSEED.  Build a list of
            # reactions, removing duplicates by checking the linked reactions field and
            # skipping obsolete reactions.
            numObsolete = 0
            unique = set()
            with open(args.rxnfilename, 'r') as handle:
                header = handle.readline()
                names = header.strip('\n').split('\t')
                if names[0] != 'id' or names[18] != 'is_obsolete' or names[19] != 'linked_reaction':
                    sys.stderr.write('Format of modelseed reaction file is invalid, column headers are wrong')
                    exit(1)
                for line in handle:
                    fields = line.strip('\n').split('\t')
                    # Skip reactions that are marked obsolete
                    if fields[18] == '1':
                        numObsolete += 1
                        continue
                    # Add the reaction and its linked reactions if they are not already
                    # in the unique set.
                    if fields[0] not in unique:
                        kbaseRxnIdList.append(fields[0])
                        unique.add(fields[0])
                        if fields[19] != 'null':
                            linkedReactions = fields[19].split(';')
                            for index in range(len(linkedReactions)):
                                unique.add(linkedReactions[index])
            print 'Found %d reactions in modelseed master reaction file' %(len(kbaseRxnIdList))
            print 'Found %d obsolete reactions in modelseed master reaction file' %(numObsolete)
            del unique
    
        else:
            print 'Source type %s is not supported' %(args.source)
            exit(1)
        
        # Find the corresponding KEGG reaction IDs for the KBase reaction IDs.  Note that not
        # all of the KBase reactions have an alias to a KEGG reaction.
        keggRxnIdList = keggExtractor.getKeggReactionList(kbaseRxnIdList, args.increment, args.biochem, args.biochemws)
        dataParser.writeKeggReactionFile(dataParser.sources['kegg']['reaction_file'], keggRxnIdList)
        print 'Stored %d aliases to KEGG reactions' %(len(keggRxnIdList))
    else:
        print 'Loading list from file "%s"...' %(dataParser.sources['kegg']['reaction_file'])
        keggRxnIdList = dataParser.readKeggReactionFile(dataParser.sources['kegg']['reaction_file'])
        print 'Retrieved %d aliases from KEGG reactions' %(len(keggRxnIdList))
    print 'Done at %s\n' %(now()) 
    step += 1
                           
    # Second step is to update the local KEGG reaction and enzyme databases if requested.
    print 'Updating local KEGG databases at %s' %(now())
    print 'Local reaction database is in "%s" and has %d records' %(args.reactionDB, keggExtractor.reactionDB.size())
    print 'Local enzyme database is in "%s" and has %d records' %(args.enzymeDB, keggExtractor.enzymeDB.size())
    if step >= args.step:
        if args.update:
            keggExtractor.updateLocalDatabases(keggRxnIdList)
        else:
            print 'KEGG database update is disabled'
    else:
        print 'KEGG database update step is skipped'
    print 'Done at %s\n' %(now())
    step += 1

    # Third step is to create a reactions to complexes file.
    print 'Creating mapping of reaction to complex at %s' %(now())
    if step >= args.step:
        print 'Saving mapping to file "%s"...' %(dataParser.sources['kegg']['reaction_complex_file'])
        reactionToComplexes = keggExtractor.mapReactionToComplex(keggRxnIdList)
        dataParser.writeReactionComplexFile(dataParser.sources['kegg']['reaction_complex_file'], reactionToComplexes)
        print 'Stored %d reaction to complex mappings' %(len(reactionToComplexes))
    else:
        print 'Loading mapping from file "%s"...' %(dataParser.sources['kegg']['reaction_complex_file'])
        reactionToComplexes = dataParser.readReactionComplexFile(dataParser.sources['kegg']['reaction_complex_file'])
        print 'Retrieved %d reaction to complex mappings' %(len(reactionToComplexes))
    print 'Done at %s\n' %(now())
    del reactionToComplexes
    step += 1

    # Fourth step is to create a complex to roles file.
    print 'Creating mapping of complex to role at %s' %(now())
    if step >= args.step:
        print 'Saving mapping to file "%s"...' %(dataParser.sources['kegg']['complex_role_file'])
        complexToRoles = keggExtractor.mapComplexToRole()
        dataParser.writeComplexRoleFile(dataParser.sources['kegg']['complex_role_file'], complexToRoles)
        print 'Stored %d complex to role mappings' %(len(complexToRoles))
    else:
        print 'Loading mappping from file "%s"' %(dataParser.sources['kegg']['complex_role_file'])
        complexToRoles = dataParser.readComplexRoleFile(dataParser.sources['kegg']['complex_role_file'])
        print 'Retrieved %d complex to role mappings' %(len(dataParser.sources['kegg']['complex_role_file']))
    print 'Done at %s\n' %(now())
    del complexToRoles
    step += 1

    # Fifth step is to create a feature ID to roles file and download amino acid sequences.
    print 'Getting amino acid sequences and feature ID to role mapping for OTU representatives at %s' %(now())
    if step >= args.step:
        print 'Saving feature ID to roles mapping to file "%s"...' %(dataParser.sources['kegg']['otu_fid_role_file'])
        print 'Saving amino acid sequences to file "%s"...' %(dataParser.sources['kegg']['protein_fasta_file'])
        fidsToRoles = keggExtractor.mapFeatureToRole(dataParser.sources['kegg']['protein_fasta_file']) #dataParser.sources['kegg']['protein_fasta_file']
        dataParser.writeFidRoleFile(dataParser.sources['kegg']['otu_fid_role_file'], fidsToRoles)
        print 'Stored %d feature ID to roles mappings' %(len(fidsToRoles))
    else:
        print 'Loading feature ID to roles mapping from file "%s"' %(dataParser.sources['kegg']['otu_fid_role_file'])
        fidsToRoles = dataParser.readFidRoleFile(dataParser.sources['kegg']['otu_fid_role_file'])
        print 'Retrieved %d feature ID to roles mappings' %(len(fidsToRoles))
    print 'Done at %s\n' %(now())
    del fidsToRoles
    step += 1

    return

if __name__ == '__main__':
    # Force output to be unbuffered.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)

    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='build-kegg-db', epilog=desc3)
    parser.add_argument('configFilePath', help='path to configuration file', action='store', default=None)
    parser.add_argument('rxnfilename', help='path to file with list of reactions', action='store')
    parser.add_argument('--source', help='type of source file', action='store', dest='source', default='template')
    parser.add_argument('--increment', help='size of increment when retrieving lists from fba modeling service', action='store', dest='increment', type=int, default=100)
    parser.add_argument('--biochem', help='ID of Biochemistry object', action='store', dest='biochem', default='default')
    parser.add_argument('--biochemws', help='ID of workspace containing Biochemistry object', action='store', dest='biochemws', default='kbase')
    parser.add_argument('--reactiondb', help='path to file with source KEGG reaction database', action='store', dest='reactionDB', default='reaction.db')
    parser.add_argument('--enzymedb', help='path to file with KEGG enzyme database', action='store', dest='enzymeDB', default='enzyme.db')
    parser.add_argument('--organismdb', help='path to file with KEGG organism database', action='store', dest='organismDB', default='organism.db')
    parser.add_argument('--no-update', help='skip adding missing reactions and enzymes to local databases', action='store_false', dest='update', default=True)
    parser.add_argument('--force', help='remove existing static database files first', dest='force', action='store_true', default=False)
    parser.add_argument('--step', help='start generating database files from this step', dest='step', action='store', type=int, default=1)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)
    print 'Generating KEGG static database files in directory "%s"' %(config['data_folder_path'])
    print 'KEGG server is at %s' %(config['kegg_url'])
    print 'FBA modeling server is at %s' %(config['fbamodeling_url'])
    print 'Biochemistry object is %s/%s' %(args.biochemws, args.biochem)
    print 'List of KBase reaction IDs comes from %s\n' %(args.rxnfilename)

    # Create a DataParser object for working with the static database files (the
    # data folder is created if it does not exist).
    dataParser = DataParser(config)

    # Create a KEGGExtractor object for getting data from KEGG.
    keggExtractor = KEGGExtractor(config, args.reactionDB, args.enzymeDB, args.organismDB)

    # Generate the intermediate static database files from the central data model.
    # Note the status file is not used when generating intermediate files.
    try:
        generate_data(dataParser, keggExtractor, args)
    except:
        sys.stderr.write('\nERROR Caught exception at %s...\n' %(now()))
        traceback.print_exc(file=sys.stderr)
        exit(1)

    exit(0) 
