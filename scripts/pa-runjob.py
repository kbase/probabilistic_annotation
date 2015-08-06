#! /usr/bin/python

import argparse
import sys
import os
import json
import traceback
from biokbase.probabilistic_annotation.Worker import ProbabilisticAnnotationWorker
from biokbase.probabilistic_annotation.Helpers import make_object_identity, timestamp, ServiceVersion, ServiceName, ProbAnnoType
from biokbase.userandjobstate.client import UserAndJobState
from biokbase.workspace.client import Workspace

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='pa-runjob')
    parser.add_argument('jobDirectory', help='path to job directory for the job', action='store', default=None)
    args = parser.parse_args()
    
    # Run the job to create a ProbAnno typed object.  Everything we need to know
    # about the job is in the jobdata.json file.
    jobDataPath = os.path.join(args.jobDirectory, 'jobdata.json')
    job = json.load(open(jobDataPath, 'r'))

    # Create a user and job state client and authenticate as the user.
    ujsClient = UserAndJobState(job['config']['userandjobstate_url'], token=job['context']['token'])
    
    try:
        # Create a Worker for running each step in the job.
        worker = ProbabilisticAnnotationWorker(job['input']['genome'], job['id'], context=job['context'])

        # Get the Genome object from the specified workspace.
        try:
            ujsClient.update_job_progress(job['id'], job['context']['token'], 'getting genome object', 1, timestamp(3600))
        except:
            pass
        wsClient = Workspace(job['config']['workspace_url'], token=job['context']['token'])
        genomeObjectId = make_object_identity(job['input']['genome_workspace'], job['input']['genome'])
        objectList = wsClient.get_objects( [ genomeObjectId ] )
        genomeObject = objectList[0]
        
        # Convert Genome object to fasta file.
        try:
            ujsClient.update_job_progress(job['id'], job['context']['token'], 'converting Genome object to fasta file', 1, timestamp(3600))
        except:
            pass
        fastaFile = worker.genomeToFasta(genomeObject['data']['features'])
        
        # Run blast using the fasta file.
        try:
            ujsClient.update_job_progress(job['id'], job['context']['token'], 'running blast', 1, timestamp(3600))
        except:
            pass
        blastResultFile = worker.runBlast(fastaFile)
        
        # Calculate roleset probabilities.
        try:
            ujsClient.update_job_progress(job['id'], job['context']['token'], 'calculating roleset probabilities', 1, timestamp(300))
        except:
            pass
        rolestringTuples = worker.rolesetProbabilitiesMarble(blastResultFile)
        
        # Build ProbAnno object and store in the specified workspace.
        try:
            ujsClient.update_job_progress(job['id'], job['context']['token'], 'building ProbAnno object', 1, timestamp(120))
        except:
            pass

        objectData = dict()
        objectData['id'] = job['input']['probanno']
        objectData['genome'] = job['input']['genome']
        objectData['genome_workspace'] = job['input']['genome_workspace'];
        objectData['roleset_probabilities'] = rolestringTuples;
        objectData['skipped_features'] = worker.skippedFeatures(genomeObject['data']['features'], blastResultFile, rolestringTuples)

        objectMetaData = dict()
        objectMetaData['num_rolesets'] = len(objectData['roleset_probabilities'])
        objectMetaData['num_skipped_features'] = len(objectData['skipped_features'])

        objectProvData = dict()
        objectProvData['time'] = timestamp(0)
        objectProvData['service'] = os.environ.get('KB_SERVICE_NAME', ServiceName)
        objectProvData['service_ver'] = ServiceVersion
        objectProvData['method'] = 'annotate'
        objectProvData['method_params'] = job['input'].items()
        objectProvData['input_ws_objects'] = [ '%s/%s/%d' %(genomeObject['info'][7], genomeObject['info'][1], genomeObject['info'][4]) ]

        objectSaveData = dict()
        objectSaveData['type'] = ProbAnnoType
        objectSaveData['name'] = job['input']['probanno']
        objectSaveData['data'] = objectData
        objectSaveData['meta'] = objectMetaData
        objectSaveData['provenance'] = [ objectProvData ]

        # Store the ProbAnno object in the specified workspace.
        retryCount = 3
        while retryCount > 0:
            try:
                objectInfo = wsClient.save_objects( { 'workspace': job['input']['probanno_workspace'], 'objects': [ objectSaveData ] } )
                break
            except HTTPError as e:
                # Hopefully this is just a temporary glitch, try again in a few seconds since we worked so hard to build the object.
                retryCount -= 1
                message = 'HTTP error %s when saving %s to workspace %s' %(e.reason, job['input']['probanno'], job['input']['probanno_workspace'])
                sys.stderr.write('[warning] '+message+'\n')
                time.sleep(15)
        
        # Saving the object failed so raise the last exception that was caught.
        if retryCount == 0:
            raise e

        # Mark the job as complete with the given status.
        ujsClient.complete_job(job['id'], job['context']['token'], 'done', None, { })

        # Cleanup resources.
        worker.cleanup()

    except Exception as e:
        # Mark the job as failed. Note that the cleanup() method is not run to leave
        # a trail for debugging the failure.
        tb = traceback.format_exc()
        sys.stderr.write(tb)
        ujsClient.complete_job(job['id'], job['context']['token'], 'failed', tb, { })
    
    exit(0)
