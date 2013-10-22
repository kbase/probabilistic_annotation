#! /usr/bin/python -u

# Note -u forces stdout to be unbuffered.

import argparse
import subprocess
import traceback
import sys
import time
from biokbase.workspaceService.Client import workspaceService
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation

''' Exception raised when a command returns a non-zero return code. '''
class CommandError(Exception):
    pass

''' Exception raised when a job ends with an error. '''
class JobError(Exception):
    pass

''' Run a command and get the output. '''

class Command:
    def __init__(self, args):
        self.proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (self.stdout, self.stderr) = self.proc.communicate()
        if self.proc.returncode != 0:
            raise CommandError(self.proc.returncode)
    
    def getLine(self, lineno):
        lines = self.stdout.split('\n')
        return lines[lineno]

''' Run the workflow for probabilistic annotation paper. '''

class Workflow:

    ''' Check if object should be created. '''

    def _isObjectMissing(self, type, id):

        # See if the object exists.
        hasObjectParams = dict()
        hasObjectParams['type'] = type
        hasObjectParams['id'] = id
        hasObjectParams['workspace'] = self.args.workspace
        hasObjectParams['auth'] = self.token
        # note - there is a has_object command, which is probably faster (since you don't have to
        # send the whole object across the pipe)
        try:
            exists = self.wsClient.get_object(hasObjectParams)
        except:
            return True
        if not exists:
            return True

        # Delete the object if it exists and force option is turned on.
        if self.args.force:
            self.wsClient.delete_object(hasObjectParams)
            self.wsClient.delete_object_permanently(hasObjectParams)
            return True

        return False

    ''' Wait for the specified job to finish. '''

    def _waitForJob(self, jobid):
        done = False
        totaltime = 0
        increment = 60
        while not done:
            time.sleep(increment)
            totaltime += 60
            jobList = self.wsClient.get_jobs( { 'jobids': [ jobid ], 'auth': self.token } )
            if jobList[0]['status'] == 'done':
                done = True
            if jobList[0]['status'] == 'error':
                print '  [ERROR]'
                print jobList[0]
                raise JobError("Job '%s' finished with an error" %(jobid))
            if totaltime > self.args.maxtime:
                print '  [ERROR]'
                print jobList[0]
                raise JobError("Job %s did not finish within specified maximum time of %d seconds" %(jobid, self.args.maxtime))
        return

    ''' Load genome into the workspace from the specified source. (Note - this step is the same for all types of Gapfill)'''

    def _loadGenome(self):
        loadGenomeParams = dict()
        loadGenomeParams['genome'] = self.args.genome
        loadGenomeParams['workspace'] = self.args.workspace
        loadGenomeParams['source'] = self.args.source
        loadGenomeParams['auth'] = self.token
        genomeMeta = self.fbaClient.genome_to_workspace(loadGenomeParams)
        return genomeMeta

    ''' Build draft model object from a genome. (Note - this step is the same for all types of Gapfill) '''

    def _buildDraftModel(self, draftModel):
        draftModelParams = dict()
        draftModelParams['genome'] = self.args.genome
        draftModelParams['genome_workspace'] = self.args.workspace
        draftModelParams['model'] = draftModel
        draftModelParams['workspace'] = self.args.workspace
        draftModelParams['auth'] = self.token
        draftModelMeta = self.fbaClient.genome_to_fbamodel(draftModelParams) 
        return draftModelMeta

    ''' Run a Probabilistic Annotation on your genome. Note - this is the same for probanno gapfill of both normal and iterative varieties ''' 

    def _runProbAnno(self, probanno):
        annotateParams = dict()
        annotateParams['genome'] = self.args.genome
        annotateParams['genome_workspace'] = self.args.workspace
        annotateParams['probanno'] = probanno
        annotateParams['probanno_workspace'] = self.args.workspace
        jobid = self.paClient.annotate(annotateParams)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        print '  Waiting for job %s to end ...' %(jobid)
        self._waitForJob(jobid)
        return

    ''' Run a Calculate (to get RxnProbs) job on your genome. Note - this is the same for probanno gapfill of both normal and iterative varieties '''
    def _runCalculate(self, rxnprobs):
        calculateParams = dict()
        calculateParams['probanno'] = self.probanno
        calculateParams['probanno_workspace'] = self.args.workspace
        calculateParams['rxnprobs'] = rxnprobs
        calculateParams['rxnprobs_workspace'] = self.args.workspace
        rxnprobsMeta = self.paClient.calculate(calculateParams)
        return rxnprobsMeta

    ''' Run gap fill on a draft model. '''

    def _gapfill(self, draftModel, model, rxnprobs, iterative=False):
        gapfillFormulation = dict()
        gapfillFormulation['directionpen'] = 4
        gapfillFormulation['singletranspen'] = 25
        gapfillFormulation['biomasstranspen'] = 25
        gapfillFormulation['transpen'] = 25
        if iterative and self.args.numsolutions > 1:
            print "  [WARNING]: Number of solutions > 1 for iterative gap fill is not allowed. We will ignore that argument for this run."
        else:
            gapfillFormulation['num_solutions'] = self.args.numsolutions

        if rxnprobs != None:
            gapfillFormulation['probabilisticReactions'] = rxnprobs
            gapfillFormulation['probabilisticAnnotation_workspace'] = self.args.workspace
        gapfillParams = dict()
        gapfillParams['model'] = draftModel
        gapfillParams['model_workspace'] = self.args.workspace
        gapfillParams['out_model'] = model
        gapfillParams['workspace'] = self.args.workspace
        gapfillParams['formulation'] = gapfillFormulation
        gapfillParams['solver'] = 'CPLEX'
        gapfillParams['auth'] = self.token
        if iterative:
            gapfillParams['completeGapfill'] = 1
            print 'turned on iterative'
        return

        job = self.fbaClient.queue_gapfill_model(gapfillParams)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        print '  Waiting for job %s to end ...' %(job['id'])
        self._waitForJob(job['id'])

    ''' Find the first solution after running gap fill. Or, if getAll is specified (needed for iterative gapfill),
        find all of the solutions after running gap fill. '''

    def _findGapfillSolution(self, model, getAll = False):
        # Get the Model object data.
        getModelParams = dict()
        getModelParams['models'] = [ model ]
        getModelParams['workspaces'] = [ self.args.workspace ]
        getModelParams['auth'] = self.token
        models = self.fbaClient.get_models(getModelParams)
        gapfillList = models[0]["unintegrated_gapfillings"]
        if len(gapfillList) < 1:
            raise IOError("Model %s/%s does not have any unintegrated gapfillings!" %(self.args.workspace, model))
        
        # Just use the first unintegrated gap fill solution.
        gapfill = gapfillList[0]
        gapfill_uuid = gapfill[1]

        # Get the hidden GapFill object attached to the Model object.
        # The solutions that we need to integrate are in there and we need to reach into 
        # it to find out how many there are and if there was a valid solution found.
        getObjectParams = dict()
        getObjectParams['id'] = gapfill_uuid
        getObjectParams['type'] = 'GapFill'
        getObjectParams['workspace'] = 'NO_WORKSPACE'
        getObjectParams['auth'] = self.token
        gapfillObject = self.wsClient.get_object(getObjectParams)
        gapfillSolutionIds = []
        for ii in range(len(gapfillObject['data']['gapfillingSolutions'])):
            gapfill_solutionid = '%s.solution.%s' %(gapfill_uuid, ii)
            gapfillSolutionIds.append(gapfill_solutionid)

        if getAll:
            return gapfillSolutionIds
        else:
            return [ gapfillSolutionIds[0] ]

    ''' Integrate the specified solutions into the specified model. solutions is an array of solutions to integrate (needed
    for iterative gapfilling) '''

    def _integrateSolutions(self, model, integratedModel, solutions, rxnprobs):
        integrateSolutionParams = dict()
        integrateSolutionParams['model'] = model
        integrateSolutionParams['model_workspace'] = self.args.workspace
        integrateSolutionParams['out_model'] = integratedModel
        integrateSolutionParams['workspace'] = self.args.workspace
        integrateSolutionParams['auth'] = self.token
        integrateSolutionParams['gapfillSolutions'] = solutions
        integrateSolutionParams['gapgenSolutions'] = [ ]
        if rxnprobs != None:
            integrateSolutionParams['rxnprobs'] = rxnprobs
            integrateSolutionParams['rxnprobs_workspace'] = self.args.workspace
        intModelMeta = self.fbaClient.integrate_reconciliation_solutions(integrateSolutionParams)

    ''' Run FBA to see if the specified model produces growth. '''

    def _runFBA(self, model, fba):
        # Note - I'm not 100% sure but I THINK you don't need to define a formulation object unless you're changing the media. Even then it might do it for you if you pass it in the parameters.
        # This needs some testing.
        runFbaParams = dict()
        runFbaParams['model'] = model
        runFbaParams['workspace'] = self.args.workspace
        runFbaParams['auth'] = self.token
        runFbaParams['fba'] = fba
        fbaMeta = self.fbaClient.runfba(runFbaParams)

        return

    def _getObjectiveValue(self, fba):
        # Get the FBA object and extract the objective.
        getObjectParams = dict()
        getObjectParams['id'] = fba
        getObjectParams['type'] = 'FBA'
        getObjectParams['workspace'] = self.args.workspace
        getObjectParams['auth'] = self.token
        fbaObject = self.wsClient.get_object(getObjectParams)
        for result in fbaObject['data']['fbaResults']:
            print '  Objective value is %s' %(result['objectiveValue'])
        return

    def __init__(self, args):
        # Save the arguments.
        self.args = args

        # Create clients for the three services.
        self.paClient = ProbabilisticAnnotation(self.args.paurl)
        self.token = self.paClient._headers['AUTHORIZATION'] # ProbAnno is an authenticated client
        self.wsClient = workspaceService(self.args.wsurl)
        self.fbaClient = fbaModelServices(self.args.fbaurl)

    ''' Run the workflow. '''

    def runStandard(self):
        print '=== Started Standard Gap Fill Workflow ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsMeta = self.wsClient.get_workspacemeta( { 'workspace': self.args.workspace, 'auth': self.token } )
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        
        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s ...' %(self.args.workspace, self.args.genome)
        print '  [OK] %s'  %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        draftModel = '%s.model' %(self.args.genome)
        print '+++ Step %d: Build draft model' %(step)
        if self._isObjectMissing('Model', draftModel):
            print '  Saving draft model to %s/%s ...' %(self.args.workspace, draftModel)
            draftModelMeta = self._buildDraftModel(draftModel)
        else:
            print '  Found draft model %s/%s' %(self.args.workspace, draftModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdModel = '%s.model.std' %(self.args.genome)
        print '+++ Step %d: Run standard gap fill on complete media' %(step)
        if self._isObjectMissing('Model', stdModel):
            print '  Submitting job and saving standard gap fill model to %s/%s ...' %(self.args.workspace, stdModel)
            self._gapfill(draftModel, stdModel, None)
        else:
            print '  Found standard gap fill model %s/%s' %(self.args.workspace, stdModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find standard gap fill unintegrated solutions (we only will integrate the first solution)' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, stdModel)
        solutionList = self._findGapfillSolution(stdModel)
        print '  Found unintegrated solution %s' %(solutionList[0])
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIntModel = "%s.model.std.int" %(self.args.genome)
        print '+++ Step %d: Integrate standard gap fill solution 0 on complete media (NOTE - you should check that the solution is optimal)' %(step)
        if self._isObjectMissing('Model', stdIntModel):
            print '  Integrating solution %s into model %s/%s ...' %(solutionList[0], self.args.workspace, stdIntModel)
            self._integrateSolutions(stdModel, stdIntModel, solutionList, None)
        else:
            print '  Found integrated standard gap fill model %s/%s' %(self.args.workspace, stdIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdCompleteFba = "%s.model.std.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of standard gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('FBA', self.stdCompleteFba):
            print '  Running fba and saving complete media standard gap fill FBA object to %s/%s' %(self.args.workspace, stdCompleteFba)
            self._runFBA(stdIntModel, stdCompleteFba)
        else:
            print '  Found complete media standard gap fill FBA object %s/%s'  %(self.args.workspace, stdCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        self._getObjectiveValue(stdCompleteFba)

        print '=== Completed Standard Gap Fill Workflow ==='

        return

    def runProbabilistic(self):
        print '=== Started Probabilistic Gap Fill Workflow ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsMeta = self.wsClient.get_workspacemeta( { 'workspace': self.args.workspace, 'auth': self.token } )
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s ...' %(self.args.workspace, self.args.genome)
        print '  [OK] %s'  %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        draftModel = '%s.model' %(self.args.genome)
        print '+++ Step %d: Build draft model' %(step)
        if self._isObjectMissing('Model', draftModel):
            print '  Saving draft model to %s/%s ...' %(self.args.workspace, draftModel)
            draftModelMeta = self._buildDraftModel(draftModel)
        else:
            print '  Found draft model %s/%s' %(self.args.workspace, draftModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probanno = '%s.probanno' %(self.args.genome)
        print '+++ Step %d: Create probabilistic annotation for genome' %(step)
        if self._isObjectMissing('ProbAnno', probanno):
            print '  Submitting job and saving probabilistic annotation to %s/%s ...' %(self.args.workspace, probanno)
            self._runProbAnno(probanno)
        else:
            print '  Found probabilistic annotation %s/%s' %(self.args.workspace, probanno)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
            
        step += 1
        rxnprobs = '%s.rxnprobs' %(self.args.genome)
        print '+++ Step %d: Calculate reaction probabilities' %(step)
        if self._isObjectMissing('RxnProbs', rxnprobs):
            print '  Saving reaction probabilities to %s/%s ...' %(self.args.workspace, rxnprobs)
            rxnprobsMeta = self._runCalculate(rxnprobs)
        else:
           print '  Found reaction probabilities %s/%s' %(self.args.workspace, rxnprobs)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probModel = '%s.model.pa' %(self.args.genome)
        print '+++ Step %d: Run probabilistic gap fill on complete media' %(step)
        if self._isObjectMissing('Model', probModel):
            print '  Submitting job and saving probabilistic gap fill model to %s/%s ...' %(self.args.workspace, probModel)
            self._gapfill(draftModel, probModel, rxnprobs)
        else:
            print '  Found probabilistic gap fill model %s/%s' %(self.args.workspace, probModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find probabilistic unintegrated solutions' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, probModel)
        solutionList = self._findGapfillSolution(probModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIntModel = "%s.model.pa.int" %(self.args.genome)
        print "+++ Step %d: Integrate probablisitic gap fill solution 0 on complete media (NOTE - you should check that the solution is optimal)" %(step)
        print '  Integrating solution %s into model %s/%s ...' %(solutionList[0], self.args.workspace, probIntModel)
        if self._isObjectMissing('Model', probIntModel):
            print '  Integrating probabilisitc gap fill solution 0 into model %s/%s ...' %(self.args.workspace, probIntModel)
            self._integrateSolutions(self.probModel, probIntModel, solutionList, rxnprobs)
        else:
            print '  Found integrated probabilistic gap fill model %s/%s' %(self.args.workspace, probIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probCompleteFba = "%s.model.pa.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of probabilistic gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('Model', probCompleteFba):
            print '  Running fba and saving complete media probabilistic gap fill FBA object to %s/%s' %(self.args.workspace, probCompleteFba)
            self._runFBA(probIntModel, probCompleteFba)
        else:
            print '  Found complete media probanno gap fill FBA object %s/%s'  %(self.args.workspace, probCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        self._getObjectiveValue(probCompleteFba)

        print '=== Completed Probabilistic Gap Fill Workflow ==='

        return

    def runIterative(self):
        print '=== Started Standard Iterative Gap Fill Workflow ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsMeta = self.wsClient.get_workspacemeta( { 'workspace': self.args.workspace, 'auth': self.token } )
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s ...' %(self.args.workspace, self.args.genome)
        print '  [OK] %s'  %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        draftModel = '%s.model' %(self.args.genome)
        print '+++ Step %d: Build draft model' %(step)
        if self._isObjectMissing('Model', draftModel):
            print '  Saving draft model to %s/%s ...' %(self.args.workspace, draftModel)
            draftModelMeta = self._buildDraftModel(self.draftModel)
        else:
            print '  Found draft model %s/%s' %(self.args.workspace, draftModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIterativeModel = '%s.model.std.iterative' %(self.args.genome)
        print '+++ Step %d: Run standard iterative gap fill on complete media' %(step)
        if self._isObjectMissing('Model', stdIterativeModel):
            print '  Submitting job and saving standard iterative gap fill model to %s/%s ...' %(self.args.workspace, stdIterativeModel)
            self._gapfill(draftModel, stdIterativeModel, None, iterative=True)
        else:
            print '  Found standard iterative gap fill model %s/%s' %(self.args.workspace, stdIterativeModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find standard iterative gap fill unintegrated solutions (we will integrate all solutions)'
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, stdIterativeModel)
        solutionList = self._findGapfillSolution(stdIterativeModel, getAll=True)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIterativeIntModel = "%s.model.std.iterative.int" %(self.args.genome)
        print '+++ Step %d: Integrate standard iterative gap fill solutions on complete media (NOTE - you should check that the solution is optimal)' %(step)
        if self._isObjectMissing('Model', stdIterativeIntModel):
            print '  Integrating iterative gapfilling solutions  into model %s/%s ...' %(self.args.workspace, stdIterativeIntModel)
            self._integrateSolutions(stdIterativeModel, stdIterativeIntModel, solutionList, None)
        else:
            print '  Found integrated standard iterative gap fill model %s/%s' %(self.args.workspace, stdIterativeIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIterativeCompleteFba = "%s.model.std.iterative.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of standard iterative gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('Model', stdIterativeCompleteFba):
            print '  Running fba and saving complete media standard iterative gap fill FBA object to %s/%s' %(self.args.workspace, stdIterativeCompleteFba)
            self._runFBA(stdIterativeIntModel, stdIterativeCompleteFba)
        else:
            print '  Found complete media standard iterative gap fill FBA object %s/%s'  %(self.args.workspace, stdIterativeCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        self._getObjectiveValue(stdIterativeCompleteFba)

        print '=== Completed Standard Iterative Gap Fill Workflow ==='

        return

    def runIterativeProbabilistic(self):
        print '=== Iterative Gap Fill Workflow (with probanno) ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsMeta = self.wsClient.get_workspacemeta( { 'workspace': self.args.workspace, 'auth': self.token } )
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s ...' %(self.args.workspace, self.args.genome)
        print '  [OK] %s'  %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        draftModel = '%s.model' %(self.args.genome)
        print '+++ Step %d: Build draft model' %(step)
        if self._isObjectMissing('Model', draftModel):
            print '  Saving draft model to %s/%s ...' %(self.args.workspace, draftModel)
            draftModelMeta = self._buildDraftModel(draftModel)
        else:
            print '  Found draft model %s/%s' %(self.args.workspace, draftModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))      

        step += 1
        probanno = '%s.probanno' %(self.args.genome)
        print '+++ Step %d: Create probabilistic annotation for genome' %(step)
        if self._isObjectMissing('ProbAnno', probanno):
            print '  Submitting job and saving probabilistic annotation to %s/%s ...' %(self.args.workspace, probanno)
            self._runProbAnno(probanno)
        else:
            print '  Found probabilistic annotation %s/%s' %(self.args.workspace, probanno)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
            
        step += 1
        rxnprobs = '%s.rxnprobs' %(self.args.genome)
        print '+++ Step %d: Calculate reaction probabilities' %(step)
        if self._isObjectMissing('RxnProbs', rxnprobs):
            print '  Saving reaction probabilities to %s/%s ...' %(self.args.workspace, rxnprobs)
            rxnprobsMeta = self._runCalculate(rxnprobs)
        else:
           print '  Found reaction probabilities %s/%s' %(self.args.workspace, rxnprobs)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        #### Above this line, everything about this is exactly the same as for probanno without iterative gapfill

        step += 1
        probIterativeModel = '%s.model.pa.iterative' %(self.args.genome)
        print '+++ Step %d: Run probabilistic iterative gap fill on complete media' %(step)
        if self._isObjectMissing('Model', probIterativeModel):
            print '  Submitting job and saving probabilistic iterative gap fill model to %s/%s ...' %(self.args.workspace, probIterativeModel)
            self._gapfill(draftModel, probIterativeModel, rxnprobs, iterative=True)
        else:
            print '  Found probabilistic gap fill model %s/%s' %(self.args.workspace, probIterativeModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find probabilistic iterative gap fill unintegrated solutions' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, probIterativeModel)
        solutionList = self._findGapfillSolution(self.probIterativeModel, getAll=True)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIterativeIntModel = "%s.model.pa.iterative.int" %(self.args.genome)
        print "+++ Step %d: Integrate probabilistic iterative gap fill solutions on complete media (NOTE - you should check that the solution is optimal)" %(step)
        print '  Integrating probabilistic iterative gap fill solutions into model %s/%s (note this takes a very long time)...' %(self.args.workspace, probIterativeIntModel)
        if self._isObjectMissing('Model', probIterativeIntModel):
            print '  Integrating probabilistic iterative gap fill solutions into model %s/%s ...' %(self.args.workspace, probIterativeIntModel)
            self._integrateSolutions(self.probIterativeModel, probIterativeIntModel, solutionList, rxnprobs)
        else:
            print '  Found integrated probabilistic iterative gap fill model %s/%s' %(self.args.workspace, probIterativeIntModell)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIterativeCompleteFba = "%s.model.pa.iterative.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of probabilistic iterative gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('Model', probIterativeCompleteFba):
            print '  Running fba and saving complete media probabilistic iterative gap fill FBA object to %s/%s' %(self.args.workspace, probIterativeCompleteFba)
            self._runFBA(self.probIterativeIntModel, probIterativeCompleteFba)
        else:
            print '  Found complete media probabilistic iterative gap fill FBA object %s/%s'  %(self.args.workspace, probIterativeCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        self._getObjectiveValue(stdIterativeCompleteFba)

        print '=== Completed Probabilistic Iterative Gap Fill Workflow ==='
        
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
    parser.add_argument('--genome_source', help='Source for genome. Valid values include "kbase", "seed", "rast". (Default is seed)', action='store', dest='source', default='seed')
    parser.add_argument('-?', '--usage', help='show usage information', action='store_true', dest='usage')
    parser.add_argument('-w', '--workspace', help='workspace for storing objects', action='store', dest='workspace')
    parser.add_argument('--force', help='Force rebuilding of all objects', action='store_true', dest='force', default=False)
    parser.add_argument('--standard', help='Run standard gap fill workflow', action='store_true', dest='standard', default=False)
    parser.add_argument('--prob', help='Run probabilistic gap fill workflow', action='store_true', dest='prob', default=False)
    parser.add_argument('--iterative', help='Run standard iterative gap fill workflow ', action='store_true', dest='iterative', default=False)
    parser.add_argument('--iterativeprob', help='Run probabilistic iterative gap fill workflow ', action='store_true', dest='iterativeprob', default=False)
    parser.add_argument('--num-solutions', help='Number of solutions to find for gap fill (ignored for iterative gap filling)', action='store', dest='numsolutions', type=int, default=10)
    parser.add_argument('--ws-url', help='URL for workspace service', action='store', dest='wsurl', default='http://www.kbase.us/services/workspace')
    parser.add_argument('--fba-url', help='URL for fba model service', action='store', dest='fbaurl', default='http://bio-data-1.mcs.anl.gov/services/fba')
    parser.add_argument('--pa-url', help='URL for probabilistic annotation service', action='store', dest='paurl', default='http://www.kbase.us/services/probabilistic_annotation')
    parser.add_argument('--knockoutdata', help='OPTIONAL. Provide a knockout data PhenotypeSet and we will create a new model gapfilled to each media in it.', action='store', default=None)
    parser.add_argument('--maxtime', 
                        help='OPTIONAL. Maximum amount of time to wait for a job to finish (by default the maximum is 2 hour for normal gapfill jobs and probanno jobs and 4 days for iterative gapfill)',
                        action='store', default=None)
    args = parser.parse_args()
    
    if args.usage:
        print parser.usage()
        exit(0)

    try:
        if args.force:
            print 'All objects will be rebuilt because --force option was specified.'
        workflow = Workflow(args)
        if args.standard:
            oldmax = args.maxtime
            if args.maxtime is None:
                args.maxtime = 7200
            workflow.runStandard()
            args.maxtime = oldmax
        if args.prob:
            oldmax = args.maxtime
            if args.maxtime is None:
                args.maxtime = 7200
            workflow.runProbabilistic()
            args.maxtime = oldmax
        if args.iterative:
            oldmax = args.maxtime
            if args.maxtime is None:
                args.maxtime = 345600
            workflow.runIterative()
            args.maxtime = oldmax
        if args.iterativeprob:
            oldmax = args.maxtime
            if args.maxtime is None:
                args.maxtime = 345600
            workflow.runIterativeProbabilistic()
            args.maxtime = oldmax
    except Exception as e:
        print '  [ERROR]'
        traceback.print_exc(file=sys.stderr)
        exit(1)

    exit(0)
        # TODO - how to check for nonzero growth?

#         if self.args.knockoutdata is not None:
#             step += 1
#             print "+++ Step %d: Gapfill to media conditions specified in the knockout data (NOT YET IMPLEMENTED)" %(step)
#             # Get a unique list of media conditions and workspaces in the phenotype object
#             # Sanity check - do we have so many media that it probably isn't a knockout phenotype dataset?
#             # Iterate over media and gapfill with normal and probanno gapfill"
# 
#             step += 1
#             print "+++ Step %d: Gapfill to media conditions specified in the knockout data (NOT YET IMPLEMENTED)" %(step)
#             # Iterate over the media again and integrate the solutions on normal and probanno gapfill (use rxnprobs for probanno)
# 
#             step += 1
#             print "+++ Step %d: Gapfill to media conditions specified in the knockout data (NOT YET IMPLEMENTED)" %(step)
#             # Iterate over the media a third time and check for growth on that media

        
