#! /usr/bin/python

import argparse
import os
import sys
import traceback
from biokbase.probabilistic_annotation.kegg.QueryKEGG import QueryKEGG
from biokbase.probabilistic_annotation.kegg.KEGGReactionDatabase import KEGGReactionDatabase
from biokbase.probabilistic_annotation.kegg.KEGGEnzymeDatabase import KEGGEnzymeDatabase
from biokbase.probabilistic_annotation.kegg.KEGGOrganismDatabase import KEGGOrganismDatabase
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
      the fba-export-template command.

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
    print 'Getting mapping of KBase reaction IDs to KEGG reaction IDs at %s' %(now())
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

    else:
        print 'Source type %s is not supported' %(args.source)
        exit(1)

    # Find the corresponding KEGG reaction IDs for the KBase reaction IDs.  Note that not
    # all of the KBase reactions have an alias to a KEGG reaction.
    keggRxnIdList = keggExtractor.getKeggReactionList(kbaseRxnIdList, args.increment, args.biochem, args.biochemws)
    dataParser.writeKeggReactionFile(dataParser.sources['kegg']['reaction_file'], keggRxnIdList)
#    keggRxnIdList = dataParser.readKeggReactionFile(dataParser.sources['kegg']['reaction_file'])
    print 'Found %d aliases to KEGG reactions' %(len(keggRxnIdList))
    print 'Done at %s\n' %(now()) 
                           
    # Second step is to update the local KEGG reaction and enzyme databases if requested.
    if args.update:
        print 'Updating local KEGG databases at %s' %(now())
        print 'Local reaction database is in "%s" and has %d records' %(args.reactionDB, keggExtractor.reactionDB.size())
        print 'Local enzyme database is in "%s" and has %d records' %(args.enzymeDB, keggExtractor.enzymeDB.size())
        keggExtractor.updateLocalDatabases(keggRxnIdList)
        print 'Done at %s\n' %(now())

    # Third step is to create a reactions to complexes file.
    print 'Creating mapping of reaction to complex at %s' %(now())
    print 'Saving reaction to complex mapping in file "%s"...' %(dataParser.sources['kegg']['reaction_complex_file'])
    reactionToComplexes = keggExtractor.mapReactionToComplex(keggRxnIdList)
    dataParser.writeReactionComplexFile(dataParser.sources['kegg']['reaction_complex_file'], reactionToComplexes)
    print 'Stored %d reaction to complex mappings' %(len(reactionToComplexes))
    print 'Done at %s\n' %(now())

    # Fourth step is to create a complex to roles file.
    print 'Creating mapping of complex to role at %s' %(now())
    print 'Saving complex to role mapping in file "%s"...' %(dataParser.sources['kegg']['complex_role_file'])
    complexToRoles = keggExtractor.mapComplexToRole()
    dataParser.writeComplexRoleFile(dataParser.sources['kegg']['complex_role_file'], complexToRoles)
    print 'Stored %d complex to role mappings' %(len(complexToRoles))
    print 'Done at %s\n' %(now())

    # Fifth step is to create a feature ID to roles file and download amino acid sequences.
    print 'Getting amino acid sequences and feature ID to role mapping for OTU representatives at %s' %(now())
    print 'Saving feature ID to roles mapping in file "%s"...' %(dataParser.sources['kegg']['otu_fid_role_file'])
    print 'Saving amino acid sequences in file "%s"...' %(dataParser.sources['kegg']['protein_fasta_file'])
    fidsToRoles = keggExtractor.mapFeatureToRole(None) #dataParser.sources['kegg']['protein_fasta_file']
    dataParser.writeFidRoleFile(dataParser.sources['kegg']['otu_fid_role_file'], fidsToRoles)
    print 'Stored %d feature ID to roles mappings' %(len(fidsToRoles))
    print 'Done at %s\n' %(now())

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
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)
    print 'Generating KEGG static database files in directory "%s"' %(config['data_folder_path'])
    print 'KEGG server is at %s' %(config['kegg_url'])
    print 'FBA modeling server is at %s' %(config['fbamodelservices_url'])
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
        sys.stderr.write("\nERROR Caught exception...\n")
        traceback.print_exc(file=sys.stderr)

    exit(0) 
