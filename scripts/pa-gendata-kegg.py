#! /usr/bin/python

import argparse
import re
import os
import sys
from biokbase.probabilistic_annotation.kegg.QueryKEGG import QueryKEGG
from biokbase.probabilistic_annotation.kegg.KEGGReactionDatabase import KEGGReactionDatabase
from biokbase.probabilistic_annotation.kegg.KEGGEnzymeDatabase import KEGGEnzymeDatabase
from biokbase.probabilistic_annotation.kegg.KEGGOrganismDatabase import KEGGOrganismDatabase
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.kegg.KEGGExtractor import KEGGExtractor
from biokbase.probabilistic_annotation.Helpers import get_config, now
from biokbase.MetabolicModelingTools.QueryDataModel import QueryDataModel, RoleNotFunctionalInFeature, NoComplexesForRole
from biokbase.fbaModelServices.Client import fbaModelServices, ServerError as FbaModelServerError

desc1 = '''
NAME
      pa-gendata-kegg -- generate static database of gene annotations from KEGG

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Generate the static database of high-quality gene annotations along with
      files containing intermediate data from KEGG.  The configFilePath argument
      specifies the path to the configuration file for the service.  The rxnfilename
      argument specifies the path to a file with a list of reactions.  The
      reactions must be defined in a Biochemistry object for the RxnProbs objects
      created by the probabilistic annotation service to be usable by the fba
      modeling service.

      The --source optional argument specifies the type of the source reaction
      file.  The supported types are "template" for a reaction file created by
      the fba-export-template command.

'''

desc3 = '''
EXAMPLES
      > pa-gendata-kegg

AUTHORS
      Mike Mundy 
'''

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
    parser.add_argument('--biochem', help='ID of Biochemistry object', action='store', dest='biochem', default=None)
    parser.add_argument('--biochemws', help='ID of workspace containing Biochemistry object', action='store', dest='biochemws', default=None)
    parser.add_argument('--reactiondb', help='path to file with source KEGG reaction database', action='store', dest='reactionDB', default='reaction.db')
    parser.add_argument('--enzymedb', help='path to file with KEGG enzyme database', action='store', dest='enzymeDB', default='enzyme.db')
    parser.add_argument('--no-update', help='skip adding missing reactions and enzymes to local databases', action='store_false', dest='update', default=True)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)

    # Create a DataParser object for working with the static database files.
    dataParser = DataParser(config)
    keggExtractor = KEGGExtractor(config, args.reactionDB, args.enzymeDB)

    print 'Generating static database files in "%s"...' %(config['data_folder_path'])
    print 'Central data model server is at %s' %(config['cdmi_url'])
    print 'FBA modeling server is at %s' %(config['fbamodelservices_url'])
    print 'KEGG server is at %s' %(config['kegg_url'])
    print 'Biochemistry object is %s/%s\n' %(args.biochemws, args.biochem)

    # First step is map KBase reaction IDs to KEGG reaction IDs.
    print 'Getting mapping of KBase reaction IDs to KEGG reaction IDs at %s' %(now())
    print 'Saving mapping to file "%s"' %(dataParser.sources['kegg']['reaction_file'])
    print 'List of KBase reaction IDs comes from %s...' %(args.rxnfilename)

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
        print 'Found %d KBase reactions in Model Template' %(len(kbaseRxnIdList))

    else:
        print '  Source type %s is not supported' %(args.source)
        exit(1)

    # Find the corresponding KEGG reaction IDs for the KBase reaction IDs.  Note that not
    # all of the KBase reactions have an alias to a KEGG reaction.
#     keggRxnIdList = keggExtractor.getKeggReactionList(kbaseRxnIdList, args.increment, args.biochem, args.biochemws)
#     dataParser.writeKeggReactionFile(dataParser.sources['kegg']['reaction_file'], keggRxnIdList)
    keggRxnIdList = dataParser.readKeggReactionFile(dataParser.sources['kegg']['reaction_file'])
    print 'Found %d aliases to KEGG reactions' %(len(keggRxnIdList))
    print 'Done at %s\n' %(now()) 
                           
    # Second step is to update the local KEGG reaction and enzyme databases if requested.
    if args.update:
        print 'Updating local KEGG databases at %s' %(now())
        print 'Local reaction database is in "%s" and has %d records...' %(args.reactionDB, keggExtractor.reactionDB.size())
        print 'Local enzyme database is in "%s" and has %d records...' %(args.enzymeDB, keggExtractor.enzymeDB.size())
        keggExtractor.updateLocalDatabases(keggRxnIdList)
        print 'Done at %s\n' %(now())

    # Third step is to create a reactions to complexes file.
    print 'Getting mapping of reaction to complexes at %s' %(now())
    print 'Saving reaction to complexes mapping in file "%s"' %(dataParser.sources['kegg']['reaction_complex_file'])
    print 'Generating from reaction and enzyme databases...'
    reactionToComplexes = keggExtractor.makeReactionComplexLinks(keggRxnIdList)
    dataParser.writeReactionComplexFile(dataParser.sources['kegg']['reaction_complex_file'], reactionToComplexes)
    print 'Stored %d reaction to complexes mappings' %(len(reactionToComplexes))
    print 'Done at %s\n' %(now())

    # Fourth step is to create a complex to roles file.
    print 'Getting mapping of complex to roles at %s' %(now())
    print 'Saving complex to roles mapping in file "%s"' %(dataParser.sources['kegg']['complex_role_file'])
    complexToRoles = keggExtractor.makeComplexRoleLinks()
    dataParser.writeComplexRoleFile(dataParser.sources['kegg']['complex_role_file'], complexToRoles)
    print 'Stored %d complex to roles mappings' %(len(complexToRoles))
    print 'Done at %s\n' %(now())

    # Fifth step is to create a feature ID to roles file and download amino acid sequences.
    print 'Getting amino acid sequences and feature ID to roles mapping for OTU representatives at %s' %(now())
    print 'Saving feature ID to roles mapping in file "%s"' %(dataParser.sources['kegg']['otu_fid_role_file'])
    print 'Saving amino acid sequences in file "%s"' %(dataParser.sources['kegg']['protein_fasta_file'])
    fidsToRoles = keggExtractor.makeFidRoleLinks(dataParser.sources['kegg']['protein_fasta_file'])
    dataParser.writeFidRoleFile(dataParser.sources['kegg']['otu_fid_role_file'], fidsToRoles)
    print 'Stored %d feature ID to roles mappings' %(len(fidsToRoles))
    print 'Done at %s\n' %(now())

    exit(0) 
