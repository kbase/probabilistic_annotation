#! /usr/bin/python

import os
import sys
import time
from biokbase.probabilistic_annotation.DataParser import getConfig, readConfig
from biokbase.auth import kb_config
from ConfigParser import ConfigParser

DefaultURL = 'https://kbase.us/services/probabilistic_annotation/'

''' Get the current URL for the service. '''

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

''' Set the current URL for the service and store in config file. '''

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
    Add deltaSeconds to the current time to get a time in the future. '''

def timestamp(deltaSeconds):
    # Just UTC timestamps to avoid timezone issues.
    now = time.time() + deltaSeconds
    ts = time.gmtime(time.time() + deltaSeconds)
    return time.strftime('%Y-%m-%dT%H:%M:%S+0000', ts)

''' Convert a job info tuple into a dictionary. '''

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
