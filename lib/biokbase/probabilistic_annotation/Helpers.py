#! /usr/bin/python

import os
import sys
import time
from biokbase.probabilistic_annotation.DataParser import getConfig, readConfig
from biokbase.auth import kb_config
from ConfigParser import ConfigParser

# Default URL for production server
DefaultURL = 'https://kbase.us/services/probabilistic_annotation/'

''' Get the current URL for the service.

    @returns Current URL string
'''

def get_url():
    # Just return the default URL when running in IRIS.
    if 'KB_RUNNING_IN_IRIS' in os.environ:
        return DefaultURL

    # Get the URL from the config file or use the default if it is not set.
    config = getConfig(kb_config)
    if 'url' in config:
        currentURL = config['url']
    else:
        currentURL = DefaultURL;
    return currentURL

''' Set the current URL for the service and store in config file.

    @param newURL New value for URL
    @returns New URL string
'''

def set_url(newURL):
    # Check for special value for the default URL.
    if newURL ==  'default':
        newURL = DefaultURL

    # Just return the URL when running in IRIS. There is no place to save it.
    if 'KB_RUNNING_IN_IRIS' in os.environ:
        return newURL

    # Save the new URL to the config file.
    config = readConfig(kb_config)
    config.set('probabilistic_annotation', 'url', newURL)
    with open(kb_config, 'w') as configfile:
        config.write(configfile)
    return newURL

''' Get a timestamp in the format required by user and job state service.

    @param deltaSeconds Seconds added to the current time to get a time in the future
    @returns Formatted timestamp string
'''

def timestamp(deltaSeconds):
    # Use UTC timestamps to avoid timezone issues.
    now = time.time() + deltaSeconds
    ts = time.gmtime(time.time() + deltaSeconds)
    return time.strftime('%Y-%m-%dT%H:%M:%S+0000', ts)

''' Convert a job info tuple into a dictionary.

    @param infoTuple Job info tuple returned by user and job state service functions
    @returns Dictionary of job info data from tuple
'''

def job_info_dict(infoTuple):
    info = dict()
    info['id'] = infoTuple[0]
    info['service'] = infoTuple[1]
    info['stage'] = infoTuple[2]
    info['started'] = infoTuple[3]
    info['status'] = infoTuple[4]
    info['last_update'] = infoTuple[5]
    info['total_progress'] = infoTuple[6]
    info['max_progress'] = infoTuple[7]
    info['progress_type'] = infoTuple[8]
    info['est_complete'] = infoTuple[9]
    info['complete'] = infoTuple[10]
    info['error'] = infoTuple[11]
    info['description'] = infoTuple[12]
    info['results'] = infoTuple[13]
    return info

''' Make an object identity structure.

    @param workspace Name or number of workspace containing object
    @param object Name or number of object
    @param ver Optional version number of object
    @returns ObjectIdentity structure for workspace APIs
'''

def make_object_identity(workspace, object, ver=None):

    objectIdentity = dict()
    if workspace.isdigit():
        objectIdentity['wsid'] = workspace
    else:
        objectIdentity['workspace'] = workspace
    if object.isdigit():
        objectIdentity['objid'] = object
    else:
        objectIdentity['name'] = object
    if ver is not None:
        objectIdentity['ver'] = ver
    return objectIdentity

''' Make working directory for a job.

    @param workDirectory Path to base working directory
    @param jobID Job identifier
    @returns Path to job directory
'''

def make_job_directory(workDirectory, jobID):
    jobDirectory = os.path.join(workDirectory, jobID)
    if not os.path.exists(jobDirectory):
        os.makedirs(jobDirectory, 0775)
    return jobDirectory
