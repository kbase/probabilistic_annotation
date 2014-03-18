#! /usr/bin/python

import os
import sys
import time
import json
import argparse
from biokbase.probabilistic_annotation.DataParser import getConfig
from biokbase.probabilistic_annotation.DataParser import DatabaseFiles
from biokbase.probabilistic_annotation.DataParser import StatusFiles
from biokbase.probabilistic_annotation.Client import _read_inifile
from biokbase.probabilistic_annotation.Shock import Client as ShockClient

desc1 = '''
NAME
      pa-savedata -- save static database of gene annotations to shock

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Save the static database of high-quality gene annotations along with
      files containing intermediate data to shock.  The files are then available
      for all servers to download.  The dataFolderPath argument specifies the
      path to the directory where the database and related files are stored.

      Note that all current instances of the file in shock are removed before
      saving the new file.  A probabilistic annotation server must be restarted
      to download and start using the new files.

      The --shock-url optional argument specifies an alternate URL for the
      Shock service.
'''

desc3 = '''
EXAMPLES
      Save static database files:
      > pa-savedata /kb/deployment/services/probabilistic_annotation/data
      
SEE ALSO
      pa-gendata

AUTHORS
      Matt Benedict, Mike Mundy 
'''

# Main script function
if __name__ == "__main__":

    # Parse arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-savedata', epilog=desc3)
    parser.add_argument('dataFolderPath', help='path to directory for storing static database files', action='store', default=None)
    parser.add_argument('--shock-url', help='url for shock service', action='store', dest='shockURL', default='http://kbase.us/services/shock-api/')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Get the user's authorization token.
    authdata = _read_inifile()
    
    # Create a shock client.
    shockClient = ShockClient(args.shockURL, authdata['token'])
    
    # Upload all of the static database files to shock.
    fileCache = { }
    for key in DatabaseFiles:
        filename = DatabaseFiles[key]
        filepath = os.path.join(args.dataFolderPath, filename)
        if os.path.exists(filepath):
            sys.stderr.write("Saving '%s'..." %(filepath))
            
            # See if the file already exists in Shock.
            query = "lookupname=ProbAnnoData/"+filename
            nodelist = shockClient.query(query)
            
            # Remove all instances of the file in Shock.
            if nodelist != None:
                for node in nodelist:
                    shockClient.delete(node["id"])
 
            # Build the attributes for this file and store as json in a separate file.
            moddate = time.ctime(os.path.getmtime(filepath))           
            attr = { "lookupname": "ProbAnnoData/"+filename, "moddate": moddate }
            attrFilename = os.path.join(args.dataFolderPath, filename+".attr")
            attrFid = open(attrFilename, "w")
            json.dump(attr, attrFid, indent=4)
            attrFid.close()
            
            # Upload the file to Shock.
            metadata = shockClient.create_node(filepath, attrFilename)
            fileCache[key] = metadata
            os.remove(attrFilename)
            
            # Remove the list of users from the read ACL to give the file public read permission.
            readacl = shockClient.get_acl(metadata["id"], "read")
            shockClient.delete_acl(metadata["id"], readacl["read"], "read")
            sys.stderr.write("done\n")
            
        else:
            sys.stderr.write("Could not find '%s' so it was not saved\n" %(filepath))
            
    # Save the metadata on all of the database files.
    cacheFilename = os.path.join(args.dataFolderPath, StatusFiles["cache_file"])
    json.dump(fileCache, open(cacheFilename, "w"), indent=4)
    
    exit(0)
