#! /usr/bin/python

import argparse
import subprocess
import traceback
import sys
import time
from biokbase.workspaceService.Client import workspaceService
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation

class CommandError(Exception):
    pass

# Make this a class.
# Method to get stdout and build array
class Command:
    def __init__(self, args):
        self.proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (self.stdout, self.stderr) = self.proc.communicate()
        if self.proc.returncode != 0:
            raise CommandError(self.proc.returncode)
    
    def getLine(self, lineno):
        lines = self.stdout.split('\n')
        return lines[lineno]

class ServiceUrl:
    def __init__(self, command, url):
        self.cmd = command
        self.url = url
        currentUrl = Command( [self.cmd] )
        self.savedUrl = currentUrl.getLine(1)
        
        newUrl = Command( [ self.cmd, self.url] )
        if newUrl.getLine(1) != self.url:
            raise CommandError("%s failed to set url to '%s'" %(self.cmd, self.url))
        
    def reset(self):
        output = run_cmd( [self.cmd, self.savedUrl] )
        lines = output.split('\n')
        if lines[1] != self.savedUrl:
            raise CommandError("%s failed to set url to '%s'" %(self.cmd, self.savedUrl))
        
def wait_for_job(jobid, wsClient, token):
    done = False
    while not done:
        time.sleep(60)
        jobList = wsClient.get_jobs( { 'jobids': [ jobid ], 'auth': token } )
        if jobList[0]['status'] == 'done':
            done = True
        if jobList[0]['status'] == 'error':
            print '  [ERROR]'
            print jobList[0]
            exit(1)        
    return
            
if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Workflow')
    parser.add_argument('genome', help='ID of genome (e.g. kb|g.422)', action='store', default=None)
    parser.add_argument('-?', '--usage', help='show usage information', action='store_true', dest='usage')
    parser.add_argument('-w', '--workspace', help='workspace for storing objects', action='store', dest='workspace')
    parser.add_argument('--num-solutions', help='number of solutions to find for gap fill', action='store', dest='numsolutions', default=10)
    parser.add_argument('--ws-url', help='url for workspace service', action='store', dest='wsurl', default='http://www.kbase.us/services/workspace')
    parser.add_argument('--fba-url', help='url for fba model service', action='store', dest='fbaurl', default='http://bio-data-1.mcs.anl.gov/services/fba')
    parser.add_argument('--pa-url', help='url for probabilistic annotation service', action='store', dest='paurl', default='http://www.kbase.us/services/probabilistic_annotation')
    args = parser.parse_args()
    
    if args.usage:
        print parser.usage()
        exit(0)

    # Create clients for the three services.
    paClient = ProbabilisticAnnotation(args.paurl)
    token = paClient._headers['AUTHORIZATION'] # ProbAnno is an authenticated client
    wsClient = workspaceService(args.wsurl)
    fbaClient = fbaModelServices(args.fbaurl)
    print token
    
    # Check up on the workspace.
    # Delete if requested.
    step = 0
    try:
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(args.paurl)
        print '  Workspace service url is %s' %(args.wsurl)
        print '  FBA model service url is %s' %(args.fbaurl)
        print '  Checking workspace %s ...' %(args.workspace)
        wsMeta = wsClient.get_workspacemeta( { 'workspace': args.workspace, 'auth': token } )
        print '  [OK]'
        
        step += 1
        print '+++ Step %d: Load genome' %(step)
        print '  Loading genome to %s/%s ...' %(args.workspace, args.genome)
        loadGenomeParams = dict()
        loadGenomeParams['genome'] = args.genome
        loadGenomeParams['workspace'] = args.workspace
        loadGenomeParams['source'] = 'kbase'
        loadGenomeParams['auth'] = token
        genomeMeta = fbaClient.genome_to_workspace(loadGenomeParams)
        print '  [OK]'
        
        step += 1
        probanno = '%s.probanno' %(args.genome)
        print '+++ Step %d: Create probabilistic annotation for genome' %(step)
        print '  Submitting job and saving probabilistic annotation to %s/%s ...' %(args.workspace, probanno)
        annotateParams = dict()
        annotateParams['genome'] = args.genome
        annotateParams['genome_workspace'] = args.workspace
        annotateParams['probanno'] = probanno
        annotateParams['probanno_workspace'] = args.workspace
        jobid = paClient.annotate(annotateParams)
        print '  [OK]'
        print '  Waiting for job %s to end ...' %(jobid)
        wait_for_job(jobid, wsClient, token)
        print '  [OK]'
        
        step += 1
        rxnprobs = '%s.rxnprobs' %(args.genome)
        print '+++ Step %d: Calculate reaction probabilities' %(step)
        print '  Saving reaction probabilities to %s/%s ...' %(args.workspace, rxnprobs)
        calculateParams = dict()
        calculateParams['probanno'] = probanno
        calculateParams['probanno_workspace'] = args.workspace
        calculateParams['rxnprobs'] = rxnprobs
        calculateParams['rxnprobs_workspace'] = args.workspace
        rxnprobsMeta = paClient.calculate(calculateParams)
        print '  [OK]'
        
        step += 1
        draftModel = '%s.model' %(args.genome)
        print '+++ Step %d: Build draft model' %(step)
        print '  Saving model to %s/%s ...' %(args.workspace, draftModel)
        draftModelParams = dict()
        draftModelParams['genome'] = args.genome
        draftModelParams['genome_workspace'] = args.workspace
        draftModelParams['model'] = draftModel
        draftModelParams['workspace'] = args.workspace
        draftModelParams['auth'] = token
        draftModelMeta = fbaClient.genome_to_fbamodel(draftModelParams)
        print '  [OK]'
        
        step += 1
        stdModel = '%s.model.std' %(args.genome)
        print '+++ Step %d: Run standard gap fill on complete media' %(step)
        print '  Submitting job and saving standard gap fill model to %s/%s ...' %(args.workspace, stdModel)
        stdGapfillFormulation = dict()
        stdGapfillFormulation['directionpen'] = 4
        stdGapfillFormulation['singletranspen'] = 25
        stdGapfillFormulation['biomasstranspen'] = 25
        stdGapfillFormulation['num_solutions'] = args.numsolutions
        stdGapfillParams = dict()
        stdGapfillParams['model'] = draftModel
        stdGapfillParams['model_workspace'] = args.workspace
        stdGapfillParams['out_model'] = stdModel
        stdGapfillParams['workspace'] = args.workspace
        stdGapfillParams['formulation'] = stdGapfillFormulation
        stdGapfillParams['solver'] = 'CPLEX'
        stdGapfillParams['auth'] = token
        job = fbaClient.queue_gapfill_model(stdGapfillParams)
        print '  [OK]'
        print '  Waiting for job %s to end ...' %(job['id'])
        wait_for_job(job['id'], wsClient, token)
        print '  [OK]'
        
        step += 1
        print '+++ Step %d: Find standard gap fill unintegrated solutions'
        print '  Checking gap fill model %s/%s ...' %(args.workspace, stdModel)
        getModelParams = dict()
        getModelParams['models'] = [ stdModel ]
        getModelParams['workspaces'] = [ args.workspace ]
        getModelParams['auth'] = token
        models = fbaClient.get_models(getModelParams)
        gapfills = models[0]["unintegrated_gapfillings"]
        if len(gapfills) < 1:
            raise IOError("Model %s/%s does not have any unintegrated gapfillings!" %(args.workspace, stdModel))
        for gapfill in gapfills:
            gapfill_uuid = gapfill[1]
            gapfill_solutionid = "%s.solution.%s" %(gapfill_uuid, options.solution)
            print "%s\t%s" %(options.modelid, gapfill_solutionid)
        
        print '  [OK]'
        
        step += 1
        probModel = '%s.model.pa' %(args.genome)
        print '+++ Step %d: Run probabilistic gap fill on complete media' %(step)
        print '  Submitting job and saving probabilistic gap fill model to %s/%s ...' %(args.workspace, probModel)
        probGapfillFormulation = dict()
        probGapfillFormulation['directionpen'] = 4
        probGapfillFormulation['singletranspen'] = 25
        probGapfillFormulation['biomasstranspen'] = 25
        probGapfillFormulation['num_solutions'] = args.numsolutions
        probGapfillFormulation['probabilisticReactions'] = rxnprobs
        probGapfillFormulation['probabilisticAnnotation_workspace'] = args.workspace
        probGapfillParams = dict()
        probGapfillParams['model'] = draftModel
        probGapfillParams['model_workspace'] = args.workspace
        probGapfillParams['out_model'] = probModel
        probGapfillParams['workspace'] = args.workspace
        probGapfillParams['formulation'] = probGapfillFormulation
        probGapfillParams['solver'] = 'CPLEX'
        probGapfillParams['auth'] = token
        job = fbaClient.queue_gapfill_model(probGapfillParams)
        print '  [OK]'
        print '  Waiting for job %s to end ...' %(job['id'])
        wait_for_job(job['id'], wsClient, token)
        print '  [OK]'
        
    except Exception as e:
        print '  [ERROR]'
        traceback.print_exc(file=sys.stderr)

    exit(0)
        