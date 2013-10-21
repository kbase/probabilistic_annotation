#! /usr/bin/python

import os
import sys
import time
import json
import argparse
from biokbase.probabilistic_annotation.DataParser import getConfig
from biokbase.probabilistic_annotation.DataParser import DatabaseFiles
from biokbase.probabilistic_annotation.DataParser import StatusFiles
from biokbase.probabilistic_annotation.Shock import Client as ShockClient

# Main script function
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='probanno-savedata')
    parser.add_argument('configFilename', help='path to config file', action='store', default=None)
    args = parser.parse_args()

    # Read the config from the file.
    config = getConfig(args.configFilename)
    
    # Get an authorization token.
    try:
        fid = open(os.path.join(os.environ["HOME"], ".kbase_auth"), "r")
        token = fid.read()
        token.strip()
        fid.close()
    except:
        sys.stderr.write("Failed to get an authorization token.\n")
        exit(1)
    
    # Create a shock client.
    shockClient = ShockClient(config["shock_url"], token)
    
    # Upload all of the static database files to shock.
    fileCache = { }
    for key in DatabaseFiles:
        filename = DatabaseFiles[key]
        filepath = os.path.join(config["data_folder_path"], filename)
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
            attrFilename = os.path.join(config["data_folder_path"], filename+".attr")
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
    cacheFilename = os.path.join(config["data_folder_path"], StatusFiles["cache_file"])
    json.dump(fileCache, open(cacheFilename, "w"), indent=4)
    
    exit(0)
