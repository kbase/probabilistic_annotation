#! /usr/bin/python

import argparse
import sys
import os
import json
import traceback
from biokbase.probabilistic_annotation.Worker import ProbabilisticAnnotationWorker
from biokbase.userandjobstate.client import UserAndJobState

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='pa-runjob')
    parser.add_argument('jobDirectory', help='path to job directory for the job', action='store', default=None)
    args = parser.parse_args()
    
    # Run the job.
    jobDataPath = os.path.join(args.jobDirectory, "jobdata.json")
    job = json.load(open(jobDataPath, 'r'))
    try:
        worker = ProbabilisticAnnotationWorker()
        worker.runAnnotate(job)
    except Exception as e:
        # Mark the job as failed.
        tb = traceback.format_exc()
        sys.stderr.write(tb)
        ujsClient = UserAndJobState(job['config']['userandjobstate_url'], token=job['context']['token'])
        ujsClient.complete_job(job['id'], job['context']['token'], 'failed', tb, { })
    
    exit(0)
