#! /usr/bin/python

import sys
import time
import os
import traceback
import argparse
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.Helpers import get_config, now
from biokbase.probabilistic_annotation.DataExtractor import *

desc1 = '''
NAME
      pa-gendata-cdm -- generate static database of gene annotations

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Generate the static database of high-quality gene annotations from the
      data files generated from the set of data sources.  The configFilePath
      argument specifies the path to the configuration file for the service.

      When the --makedb optional argument is specified only build the search
      database for the configured search program.  Note that the input protein
      FASTA file must be available before using this option.

      When the --force optional argument is specified all existing files are
      deleted before they are generated.
'''

desc3 = '''
EXAMPLES
      Generate static database files:
      > pa-gendata gendata.cfg
      
SEE ALSO
      pa-gendata-cdm
      pa-gendata-kegg
      pa-savedata
      pa-loaddata

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

def generate_data(dataParser, config, force):
    
    # When regenerating the database files, remove all of them first.
    if force:
        print('Removing all static database files...')
        for filename in dataParser.DataFiles.values():
            safeRemove(filename)
        print('Done at %s' %(now()))
    
    print('Generating static database files in "%s"...\n' %(config['data_folder_path']))
    
    # Merge the data files from all of the configured sources.
    print('Merging static database files from configured sources at %s' %(now()))
    print('Sources are %s ...' %(config['data_sources']))
    dataParser.mergeDataFiles()
    print('Done at %s\n' %(now()))

    # Make the database used by the search program.
    print('Making search database for search program at %s' %(now()))
    print('Search program is %s ...' %(config['search_program']))
    dataParser.buildSearchDatabase()
    print('Done at %s\n' %(now()))
    
    print('Done generating static database files')
    return

# Main script function
if __name__ == "__main__":
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

    # Update the status file.
    dataParser.writeStatusFile('building')

    # Generate the static database files.
    try:
        if args.makeDB:
            dataParser.buildSearchDatabase()
        else:
            generate_data(dataParser, config, args.force)
        status = "ready"
    except:
        status = "failed"
        sys.stderr.write("\nERROR Caught exception...\n")
        traceback.print_exc(file=sys.stderr)
    
    # Update the status file.
    dataParser.writeStatusFile(status)
    
    exit(0)
