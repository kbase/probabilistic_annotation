#! /usr/bin/python

import argparse
import os
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.probabilistic_annotation.Helpers import get_config
from biokbase import log

desc1 = '''
NAME
      pa-loaddata -- load static database of gene annotations

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Load the static database of high-quality gene annotations along with
      files containing intermediate data.  The files are then available for
      a probabilistic annotation server on this system.  Since downloading
      from Shock can take a long time, run this command to load the static
      database files before the server is started.  The configFilePath argument
      specifies the path to the configuration file for the service.

      Note that a probabilistic annotation server is unable to service client
      requests for the annotate() and calculate() methods while this command is
      running and must be restarted to use the new files.
'''

desc3 = '''
EXAMPLES
      Load static database files:
      > pa-loaddata loaddata.cfg
      
SEE ALSO
      pa-gendata
      pa-savedata

AUTHORS
      Matt Benedict, Mike Mundy 
'''

# Main script function
if __name__ == "__main__":

    # Parse arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-loaddata', epilog=desc3)
    parser.add_argument('configFilePath', help='path to configuration file', action='store', default=None)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create a log object.
    submod = os.environ.get('KB_SERVICE_NAME', 'probabilistic_annotation')
    mylog = log.log(submod, ip_address=True, authuser=True, module=True, method=True,
        call_id=True, config=args.configFilePath)

    # Get the probabilistic_annotation section from the configuration file.
    config = get_config(args.configFilePath)

    # Create a DataParser object for working with the static database files (the
    # data folder is created if it does not exist).
    dataParser = DataParser(config)

    # Get the static database files.  If the files do not exist and they are downloaded
    # from Shock, the command may run for a long time.
    testDataPath = os.path.join(os.environ['KB_TOP'], 'services', submod, 'testdata')
    dataOption = dataParser.getDatabaseFiles(mylog, testDataPath)

    exit(0)
