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
    description=''' Workflow.py

    DESCRIPTION:
           Performs the following steps in our analysis pipeline:
           1: Load a genome from the specified source (default: SEED) into the workspace
           2: Generate probabilistic annotation for the genome
           3: Calcualte rxnprobs from the probabilistic annotation
           4: Generate a draft model (ungapfilled) from the genome
           5: Run non-probanno gapfill (on complete media)
           6: Run probanno gapfill (on complete media)
           7: Integrate the non-probanno gapfill solution on complete media (no GPRs)
           8: Integrate the probanno gapfill on complete media (including GPRs)
           9: Check for non-zero growth for both probanno and non-probannog gapfill. Report failure if growth < 1E-5
           IF a knockout phenotype is provided:
           FOR EACH knockout phenotype media condition:
               9: Gapfill to that media
              10: Integrate the solution
              11: Check for growth on that media. Report failure if growth < 1E-5.
    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Workflow', description=description)
    parser.add_argument('genome', help='ID of genome that matches the source (e.g. 83333.1 for SEED, kb|g.0 for kbase)', action='store', default=None)
    parser.add_argument('-?', '--usage', help='show usage information', action='store_true', dest='usage')
    parser.add_argument('-w', '--workspace', help='workspace for storing objects', action='store', dest='workspace')
    parser.add_argument('--genome_source', help='Source for genome. Valid values include "kbase", "seed", "rast". (Default is seed)', action='store', default='seed')
    parser.add_argument('--num-solutions', help='Number of solutions to find for gap fill', action='store', dest='numsolutions', default=10)
    parser.add_argument('--ws-url', help='URL for workspace service', action='store', dest='wsurl', default='http://www.kbase.us/services/workspace')
    parser.add_argument('--fba-url', help='URL for fba model service', action='store', dest='fbaurl', default='http://bio-data-1.mcs.anl.gov/services/fba')
    parser.add_argument('--pa-url', help='URL for probabilistic annotation service', action='store', dest='paurl', default='http://www.kbase.us/services/probabilistic_annotation')
    parser.add_argument('--knockoutdata', help='OPTIONAL. Provide a knockout data PhenotypeSet and we will create a new model gapfilled to each media in it.', action='store', default=None)
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
        print '+++ Step %d: Load genome from the specified source' %(step)
        print '  Loading genome to %s/%s ...' %(args.workspace, args.genome)
        loadGenomeParams = dict()
        loadGenomeParams['genome'] = args.genome
        loadGenomeParams['workspace'] = args.workspace
        loadGenomeParams['source'] = args.source
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
        stdGapfillFormulation['transpen'] = 25
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
        # This part I think isn't quite right. But it might be a start. Should look in my script to get gapfill solution IDs, 
        # I think there are some caveats like sometimes the first gapfill solution is empty for whatever reason!
#        for gapfill in gapfills:
#            gapfill_uuid = gapfill[1]
#            gapfill_solutionid = "%s.solution.%s" %(gapfill_uuid, options.solution)
#            print "%s\t%s" %(options.modelid, gapfill_solutionid)
        # This is a placeholder so I know the name of the variable that will be saved.
        stdSolutionToIntegrate = None
        
        print '  [OK]'
        
        step += 1
        probModel = '%s.model.pa' %(args.genome)
        print '+++ Step %d: Run probabilistic gap fill on complete media' %(step)
        print '  Submitting job and saving probabilistic gap fill model to %s/%s ...' %(args.workspace, probModel)
        probGapfillFormulation = dict()
        probGapfillFormulation['directionpen'] = 4
        probGapfillFormulation['singletranspen'] = 25
        probGapfillFormulation['biomasstranspen'] = 25
        probGapfillFormulation['transpen'] = 25
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

        # Will need to generate this one too. Procedure should be the same and we could probably use a function to do this in general.
        probSolutionToIntegrate = None

        step += 1
        stdIntModel = "%s.model.std.int" %(args.genome)
        print "+++ Step %d: Integrate standard gapfilling solution 0 on complete media (NOTE - you should check that the solution is optimal)" %(step)
        stdIntegrateSolutionParams = dict()
        stdIntegrateSolutionParams['model'] = stdModel
        stdIntegrateSolutionParams['model_workspace'] = workspace
        stdIntegrateSolutionParams['out_model'] = stdIntModel
        stdIntegrateSolutionParams['workspace'] = workspace
        stdIntegrateSolutionParams['auth'] = token
        stdIntegrateSolutionParams['gapfillSolutions'] = [ stdSolutionToIntegrate ]
        stdIntModelMeta = fbaClient.integrate_reconciliation_solutions(stdIntegrateSolutionParams)

        step += 1
        probIntModel = "%s.model.pa.int" %(args.genome)
        print "+++ Step %d: Integrate probanno gapfilling solution 0 on complete media (NOTE - you should check that the solution is optimal)" %(step)
        probIntegrateSolutionParams = dict()
        probIntegrateSolutionParams['model'] = probModel
        probIntegrateSolutionParams['model_workspace'] = workspace
        probIntegrateSolutionParams['out_model'] = probIntModel
        probIntegrateSolutionParams['workspace'] = workspace
        probIntegrateSolutionParams['auth'] = token
        probIntegrateSolutionParams['gapfillSolutions'] = [ probSolutionToIntegrate ]
        # Important! Include the rxnprobs! If this fails, it means we're using an old version of fbaModelServices. Which is bad.
        probIntegrateSolutionParams['rxnprobs'] = rxnprobs
        probIntegrateSolutionParams['rxnprobs_workspace'] = workspace
        probIntModelMeta = fbaClient.integrate_reconciliation_solutions(probIntegrateSolutionParams)

        step += 1
        stdCompleteFba = "%s.model.std.int.fba" %(args.genome)
        print "+++ Step %d: a: Check for growth of standard gapfill model on complete media (all available transporters to the cell are turned on)" %(step)
        # Note - I'm not 100% sure but I THINK you don't need to define a formulation object unless you're changing the media. Even then it might do it for you if you pass it in the parameters.
        # This needs some testing.
        stdCompleteFbaParams = dict()
        stdCompleteFbaParams['model'] = stdIntModel
        stdCompleteFbaParams['workspace'] = workspace
        stdCompleteFbaParams['auth'] = token
        stdFbaMeta = fbaClient.runfba(stdCompleteFbaParams)
        # TODO - how to check nonzero growth? I THINK this is somewhere in the metadata.

        probCompleteFba = "%s.model.pa.int.fba" %(args.genome)
        print "+++ Step %d: b: Check for growth of probanno gapfill model on complete media (all available transporters to the cell are turned on)" %(step)
        probCompleteFbaParams = dict()
        probCompleteFbaParams['model'] = probIntModel
        probCompleteFbaParams['workspace'] = workspace
        probCompleteFbaParams['auth'] = token
        probFbaMeta = fbaClient.runfba(probCompleteFbaParams)
        # TODO - how to check for nonzero growth?

        if args.knockoutdata is not None:
            step += 1
            print "+++ Step %d: Gapfill to media conditions specified in the knockout data (NOT YET IMPLEMENTED)" %(step)
            # Get a unique list of media conditions and workspaces in the phenotype object
            # Sanity check - do we have so many media that it probably isn't a knockout phenotype dataset?
            # Iterate over media and gapfill with normal and probanno gapfill"

            step += 1
            print "+++ Step %d: Gapfill to media conditions specified in the knockout data (NOT YET IMPLEMENTED)" %(step)
            # Iterate over the media again and integrate the solutions on normal and probanno gapfill (use rxnprobs for probanno)

            step += 1
            print "+++ Step %d: Gapfill to media conditions specified in the knockout data (NOT YET IMPLEMENTED)" %(step)
            # Iterate over the media a third time and check for growth on that media

    except Exception as e:
        print '  [ERROR]'
        traceback.print_exc(file=sys.stderr)

    exit(0)
        
