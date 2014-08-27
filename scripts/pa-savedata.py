#! /usr/bin/python

import argparse
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.Helpers import get_config

desc1 = '''
NAME
      pa-savedata -- save static database of gene annotations to shock

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Save the static database of high-quality gene annotations along with
      files containing intermediate data to Shock.  The files are then available
      for all servers to download.  The configFilePath argument specifies the
      path to the configuration file for the service.

      Note that all current instances of the file in Shock are removed before
      saving the new file.  A probabilistic annotation server must be restarted
      to download and start using the new files.
'''

desc3 = '''
EXAMPLES
      Save static database files:
      > pa-savedata gendata.cfg
      
SEE ALSO
      pa-gendata

AUTHORS
      Matt Benedict, Mike Mundy 
'''

# Main script function
if __name__ == "__main__":

    # Parse arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-savedata', epilog=desc3)
    parser.add_argument('configFilePath', help='path to configuration file', action='store', default=None)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)

    # Create a DataParser object for working with the static database files.
    dataParser = DataParser(config)

    # Store the static database files in Shock.
    dataParser.storeDatabaseFiles()

    exit(0)
