#! /usr/bin/python -u

# Note -u forces stdout to be unbuffered.

import argparse
import getpass
import re
import subprocess
import traceback
import sys
import time
from operator import itemgetter
from biokbase.workspace.client import Workspace
from biokbase.workspaceService.client import workspaceService as oldWorkspace
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation, _read_inifile
from biokbase.probabilistic_annotation.Helpers import job_info_dict
from biokbase.userandjobstate.client import UserAndJobState, ServerError as JobStateServerError

''' Exception raised when a command returns a non-zero return code. '''
class CommandError(Exception):
    pass

''' Exception raised when a job ends with an error. '''
class JobError(Exception):
    pass

''' Exception raised when required data is not found when trying to integrate phenotype data into already-existing gapfill solutions '''
class NoDataError(Exception):
    pass

''' Exception raised when something goes wrong with gapfill and it doesn't give us growth when we need it to '''
class NoGrowthError(Exception):
    pass

''' Exception raised when we cannot find a workspace '''
class CannotFindWorkspaceError(Exception):
    pass

''' Exception raised when unable to get an object from a workspace '''
class CannotGetObjectError(Exception):
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

    def _checkIfWorkspaceExists(self, workspace):
        workspaceIdentity = dict()
        workspaceIdentity["workspace"] = workspace
        try:
            self.wsClient.get_workspace_info( workspaceIdentity )
        except:
            raise CannotFindWorkspaceError("Could not find workspace %s. Do you have permission to read and write to it?")

    def _makeObjectIdentity(self, type, id, workspace=None, idtype = 'name'):
        # Create a valid ObjectIdentity object
        objectIdentity = dict()

        # This function accepts multiple input formats but we for now only
        # need to support two: name/workspace pairs and references.
        if idtype == 'ref':
            # The get_models function adds extra stuff to the beginning of the ID that causes get_objects to choke.
            # This regex gets rid fo that.
            good_ref = re.search(r'\d+/\d+/\d+', id)
            assert (good_ref is not None)
            objectIdentity['ref'] = good_ref.group(0)
        elif idtype == 'name':
            if workspace is None:
                objectIdentity['workspace'] = self.args.workspace
            else:
                objectIdentity['workspace'] = workspace
            objectIdentity['name'] = id
        return objectIdentity


    ''' Get an object of a certain type with a certain ID '''
    def _getWsObject(self, type, id, workspace=None, idtype='name'):      
        objectIdentity = self._makeObjectIdentity(type, id, workspace=workspace, idtype = idtype)
        try:
            objs = self.wsClient.get_objects( [ objectIdentity ] )
        except:
            raise CannotGetObjectError("Unable to get %s object %s/%s" %(type, self.args.workspace, id) )
        return objs[0]        

    ''' Check if object exists. If so, get the objectinfo, which includes name, ref, id, etc.
    Note this is the replacement for hasobject() in the old workspace server.'''
    def _getObjectInfo(self, type, id, workspace=None, idtype = 'name'):
        objectIdentity = self._makeObjectIdentity(type, id, workspace=workspace, idtype=idtype)
        # The get_object_info() syntax is set to change in a few months. Lets try to avoid this failing
        # because of a trivial syntax change...
        try:
            objinfo = self.wsClient.get_object_info( [ objectIdentity ], 0 )
        except:
            # Note that ServerError could mean either bad syntax (i.e. they replaced the function with another of the same name)
            # or the object doesn't exist. Great...
            get_obj_info_params = dict()
            get_obj_info_params['objects'] = [ objectIdentity ]
            get_obj_info_params['includeMetadata'] = 0
            get_obj_info_params['ignoreErrors'] = 0
            try:
                objinfo = self.wsClient.get_object_info(get_obj_info_params)
            except:
                try:
                    objinfo = self.wsClient.get_object_info_new(get_obj_info_params)
                except:
                    raise CannotGetObjectError("Unable to get %s object %s/%s" %(type, self.args.workspace, id) )

        return objinfo[0]

    ''' Check if object should be created. '''
    def _isObjectMissing(self, type, id, workspace=None):

        try:
#            exists = self._getWsObject(type, id, workspace=workspace, idtype='name')
            exists = self._getObjectInfo(type, id, workspace=workspace, idtype='name')
        except CannotGetObjectError:
            return True
        if not exists:
            return True

        # Delete the object if it exists and force option is turned on.
        if self.args.force:
            identity = self._makeObjectIdentity(type, id, workspace=workspace, idtype='name')
            self.wsClient.delete_objects( [ identity ] )
            return True

        return False

    ''' Wait for the specified job to finish. '''
    # Note - this is DEPRECIATED but we still keep it because it's needed for FBA.
    # Eventually everything will switch to UJS and we'll use waitForUjs for all jobs.
    def _waitForJob(self, jobid):
        done = False
        totaltime = 0
        increment = 60
        while not done:
            time.sleep(increment)
            totaltime += increment
            jobList = self.jobClient.get_jobs( { 'jobids': [ jobid ], 'auth': self.token } )
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

    '''Wait for a job submitted using the User and Job State service'''
    # Note - this should be used for probanno jobs. FBA jobs still use the old queue management system for now.
    def _waitForUjs(self, jobid):
        done = False
        totaltime = 0
        increment = 60
        while not done:
            time.sleep(increment)
            totaltime += increment

            jobinfo = None
            try:
                jobinfo = self.ujsClient.get_job_info(jobid)
            except:
                print '  [ERROR]'
                raise JobError('Unable to check status of job %s' %(jobid))

            infodict = job_info_dict(jobinfo)

            # Check if the job had an error.
            if infodict['error']:
                print '  [ERROR]'
                print "Job '%s' (%s) ended with error '%s'." %(infodict['id'], infodict['description'], infodict['status'])
                print 'Error details:'
                print self.ujsClient.get_detailed_error(infodict['id'])
                self.ujsClient.delete_job(infodict['id'])
                raise JobError('Error')
            # Check if the job is complete.
            elif infodict['complete']:
                print "Job '%s' (%s) completed successfully." %(infodict['id'], infodict['description'])
                self.ujsClient.delete_job(infodict['id'])
                done = True
                break
            # Job is still running.
            else:
                print "Job '%s' (%s) has status '%s' and is working on task %s of %s." \
                    %(infodict['id'], infodict['description'], infodict['status'], infodict['total_progress'], infodict['max_progress'])

            if totaltime > self.args.maxtime:
                print '  [ERROR]'
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
        # Probanno annotate jobs use the user and job state service.
        self._waitForUjs(jobid)
        return

    ''' Run a Calculate (to get RxnProbs) job on your genome. Note - this is the same for probanno gapfill of both normal and iterative varieties '''
    def _runCalculate(self, probanno, rxnprobs):
        calculateParams = dict()
        calculateParams['probanno'] = probanno
        calculateParams['probanno_workspace'] = self.args.workspace
        calculateParams['rxnprobs'] = rxnprobs
        calculateParams['rxnprobs_workspace'] = self.args.workspace
        rxnprobsMeta = self.paClient.calculate(calculateParams)
        return rxnprobsMeta

    ''' Run gap fill on a draft model. '''
    def _gapfill(self, draftModel, model, rxnprobs, iterative=False, media=None, mediaws="KBaseMedia"):
        gapfillFormulation = dict()

        # From looking at 4 organisms, you get the best accuracy by using a fairly high directionality penalty (8 or 12 as opposed to 4)
        # and a low to moderate transporter penalty (15 or 25). I use these since they also give good consistency with B subtilis knockout tests.
        gapfillFormulation['directionpen'] = 12
        gapfillFormulation['singletranspen'] = 25
        gapfillFormulation['biomasstranspen'] = 25
        gapfillFormulation['transpen'] = 25

        gapfillFormulation['allowedcmps'] = [ 'c', 'e', 'p' ]
        gapfillFormulation['nobiomasshyp'] = 1
        gapfillFormulation['formulation'] = {}

        if media is not None:
            gapfillFormulation['formulation']['media'] = media
            gapfillFormulation['formulation']['media_workspace'] = mediaws

        if iterative:
            # This is necessary to prevent ModelSEED from defaulting to just optimizing biomass.
            gapfillFormulation['formulation']['objectiveTerms'] = []

        if iterative and self.args.numsolutions > 1:
            print "  [WARNING] Number of solutions > 1 for iterative gap fill is not allowed. We will ignore that argument for this run."
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
#        gapfillParams['totalTimeLimit'] = self.args.maxtime
#        gapfillParams['timePerSolution'] = int(self.args.maxtime/self.args.numsolutions)
        if iterative:
            gapfillParams['completeGapfill'] = '1'
            # This is necessary until we get the automatic integration pushed to the SEED servers
            gapfillParams['integrate_solution'] = '1'
#            gapfillParams['timePerSolution'] = self.args.maxtime

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
        gapfill_id = gapfill[0]

        # Get the hidden GapFill object attached to the Model object.
        # The solutions that we need to integrate are in there and we need to reach into 
        # it to find out how many there are and if there was a valid solution found.
        gapfillObject = self._getWsObject('Gapfill', gapfill_id, workspace=self.args.workspace, idtype='name')

        gapfillSolutionIds = []
        for ii in range(len(gapfillObject['data']['gapfillingSolutions'])):
            # Recent change: Gf solution IDs start counting from 1 now
            gapfill_solutionid = '%s.gfsol.%s' %(gapfill_id, ii + 1)
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

    def _runFBA(self, model, fba, media=None, mediaws="KBaseMedia"):
        # Note - I'm not 100% sure but I THINK you don't need to define a formulation object unless you're changing the media.
        # This needs some testing.
        #
        # If no media is provided the default is complete media (everything with transporters is allowed into the cell)
        runFbaParams = dict()
        runFbaParams['model'] = model
        runFbaParams['workspace'] = self.args.workspace
        runFbaParams['auth'] = self.token
        runFbaParams['fba'] = fba
        if media is not None:
            runFbaParams['formulation'] = {}
            runFbaParams['formulation']['media'] = media
            runFbaParams['formulation']['media_workspace'] = mediaws

        fbaMeta = self.fbaClient.runfba(runFbaParams)

        return fbaMeta

    def _getObjectiveValue(self, fba):
        # Get the FBA object and extract the objective.
        # (note - the objective value is in the metadata so we should be getting this from _runFBA to not run into troubles with multiple FBAs in a single model)
        fbaObject = self._getWsObject('FBA', fba)
        # maybe should save all of these. But I don't know what to do with them.
        value = fbaObject['data']['objectiveValue']
        print '  Objective value is %s' %(value)
        return float(value)

    def _getMediaIdsFromPhenotypeSet(self, phenoid, workspace=None):
        # Get a list of media IDs and workspaces from a phenotypeSet object
        # Returns a list of unique (media, workspace) tuples
        phenotypeSetObject = self._getWsObject('PhenotypeSet', phenoid, workspace=workspace)

        mediaRefs = set()
        for phenotypeSet in phenotypeSetObject['data']['phenotypes']:
            for phenotype in phenotypeSet:
                mediaRefs.add( (phenotypeSet["media_ref"] ) )

        # Convert references to (workspace, ID) tuples so that we can call gap fill/FBA with it.
        mediaTuples = set()
        for ref in mediaRefs:
            objinfo = self._getObjectInfo('Media', ref, workspace=None, idtype = 'ref')
            mediaTuples.add( (objinfo[1], objinfo[7]) )

        return mediaTuples

    ''' Simulate a PhenotypeSet against a model, resulting in a simulation set. Note that it is quite possible that the
    phenotypeset is not in self.args.workspace. There's a sort-of standard place for this but not as much as for media.'''

    def _simulatePhenotype(self, model, phenoset, phenows, simset, positive_transporters = 0, all_transporters = 0):
        simulatePhenotypeParams = dict()
        simulatePhenotypeParams['model'] = model
        simulatePhenotypeParams['model_workspace'] = self.args.workspace
        simulatePhenotypeParams['phenotypeSet'] = phenoset
        simulatePhenotypeParams['phenotypeSet_workspace'] = phenows
        simulatePhenotypeParams['phenotypeSimultationSet'] = simset
        simulatePhenotypeParams['workspace'] = self.args.workspace
        simulatePhenotypeParams['auth'] = self.token

        simulatePhenotypeParams['positive_transporters'] = positive_transporters
        simulatePhenotypeParams['all_transporters'] = all_transporters
            
        objmeta = self.fbaClient.simulate_phenotypes(simulatePhenotypeParams)
        
        return objmeta

    ''' Find the ID for an integrated iterative gapfill to pass to the reaction sensntivity analysis. '''
    def _findIterativeGapfillUUID(self, model):       
        getModelParams = dict()
        getModelParams['models'] = [ model ]
        getModelParams['workspaces'] = [ self.args.workspace ]
        getModelParams['auth'] = self.token
        models = self.fbaClient.get_models(getModelParams)
        gapfillList = models[0]["integrated_gapfillings"]

        for gfarray in gapfillList:
            gapfill = gfarray[0]
            gapfill_obj = self._getWsObject('Gapfill', gapfill, workspace=self.args.workspace, idtype='name')
            if int(gapfill_obj["data"]["completeGapfill"]) == 1:
                # I am assuming here that there is only one iterative gapfill done on the model. When using
                # this workflow, that is a valid assumption.
                return gapfill
        return None
    
    ''' Run a reaction sensntivity analysis, saving the results in the specified rxnsensitivity object.
        IMPORTANT: Reactions are tested (deleted) in the order in which they are given in reactions_to_delete.'''
    def _runReactionSensitivity(self, model, rxnsensitivity, rxnprobs=None, deleteReactions=True):
        gapfillUUID = self._findIterativeGapfillUUID(model)
        if gapfillUUID is None:
            raise IOError("For some reason we couldnt find an integrated iterative gapfill UUID")

        reactionSensitivityParams = dict()
        reactionSensitivityParams['model'] = model
        reactionSensitivityParams['model_ws'] = self.args.workspace
        reactionSensitivityParams['rxnsens_uid'] = rxnsensitivity
        reactionSensitivityParams['workspace'] = self.args.workspace
        reactionSensitivityParams['auth'] = self.token
        reactionSensitivityParams['gapfill_solution_id'] = "%s.gfsol.1" %(gapfillUUID)

        if rxnprobs is not None:
            reactionSensitivityParams['rxnprobs_id'] = rxnprobs
            reactionSensitivityParams['rxnprobs_ws'] = self.args.workspace

        if deleteReactions:
            reactionSensitivityParams['delete_noncontributing_reactions'] = 1
        else:
            reactionSensitivityParams['delete_noncontributing_reactions'] = 0
        
        # Reaction sensitivity is a long-running job.
        # Since it does not queue by default we need to do it manually to avoid timeouts.
        queue_job_params = dict()
        queue_job_params['method'] = 'reaction_sensitivity_analysis'
        queue_job_params['parameters'] = reactionSensitivityParams
        job = self.fbaClient.queue_job(queue_job_params)
        
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        print '  Waiting for job %s to end ...' %(job['id'])
        self._waitForJob(job['id'])

        return

    ''' Run a delete_noncontributing_reactions job to delete unnecessary gapfill reactions as identified by
    a reaction sensitivity analysis (only relevant for iterative gapfill) '''
    def _runReactionDeletion(self, model, filteredmodel, rxnsensitivity):
        deleteReactionsParams = dict()
        deleteReactionsParams['new_model_uid'] = filteredmodel
        deleteReactionsParams['workspace'] = self.args.workspace
        deleteReactionsParams['rxn_sensitivity'] = rxnsensitivity
        deleteReactionsParams['rxn_sensitivity_ws'] = self.args.workspace
        deleteReactionsParams['auth'] = self.token

        objmeta = self.fbaClient.delete_noncontributing_reactions(deleteReactionsParams)

        return objmeta

    def __init__(self, args):
        # Save the arguments.
        self.args = args

        # Create clients for the three services and the jobs [temporary]
        self.paClient = ProbabilisticAnnotation(self.args.paurl)
        self.token = self.paClient._headers['AUTHORIZATION'] # ProbAnno is an authenticated client
        self.wsClient = Workspace(self.args.wsurl, token=self.token)
        self.jobClient = oldWorkspace(self.args.joburl, token=self.token)
        self.fbaClient = fbaModelServices(self.args.fbaurl, token=self.token)
        self.ujsClient = UserAndJobState(self.args.ujsurl, token=self.token)

    ''' Run the workflow. '''

    def runStandard(self):
        print '=== Started Network-based Gap Fill Workflow ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsMeta = self._checkIfWorkspaceExists(self.args.workspace)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        
        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s' %(self.args.workspace, self.args.genome)
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
        print '+++ Step %d: Run network-based gap fill on complete media' %(step)
        if self._isObjectMissing('Model', stdModel):
            print '  Submitting job and saving network-based gap fill model to %s/%s ...' %(self.args.workspace, stdModel)
            self._gapfill(draftModel, stdModel, None)
        else:
            print '  Found standard gap fill model %s/%s' %(self.args.workspace, stdModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find network-based gap fill unintegrated solutions (we only will integrate the first solution)' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, stdModel)
        solutionList = self._findGapfillSolution(stdModel)
        print '  Found unintegrated solution %s' %(solutionList[0])
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIntModel = "%s.model.std.int" %(self.args.genome)
        print '+++ Step %d: Integrate network-based gap fill solution 1 on complete media (NOTE - you should check that the solution is optimal)' %(step)
        if self._isObjectMissing('Model', stdIntModel):
            print '  Integrating solution %s into model %s/%s ...' %(solutionList[0], self.args.workspace, stdIntModel)
            self._integrateSolutions(stdModel, stdIntModel, solutionList, None)
        else:
            print '  Found integrated network-based gap fill model %s/%s' %(self.args.workspace, stdIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdCompleteFba = "%s.model.std.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of network-based gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('FBA', stdCompleteFba):
            print '  Running fba and saving complete media standard gap fill FBA object to %s/%s' %(self.args.workspace, stdCompleteFba)
            self._runFBA(stdIntModel, stdCompleteFba)
        else:
            print '  Found complete media network-based gap fill FBA object %s/%s'  %(self.args.workspace, stdCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        obj = self._getObjectiveValue(stdCompleteFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Network-based gapfilled model did not grow on minimal media after gapfilling'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        step += 1
        stdMinimalModel = "%s.model.std.int.minimal" %(self.args.genome)
        print "+++ Step %d: Gapfill to minimal media (Carbon-D-Glucose)" %(step)
        if self._isObjectMissing('Model', stdMinimalModel):
            print '   Gapfilling to Minimal media (Carbon-D-Glucose)'
            self._gapfill(stdIntModel, stdMinimalModel, None, media='Carbon-D-Glucose')
        else:
            print '   Found carbon-D-glucose gapfilled model %s/%s' %(self.args.workspace, stdMinimalModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdMinimalIntModel = "%s.model.std.int.minimal.int" %(self.args.genome)
        print "+++ Step %d: Integrate gapfill solution to minimal media (Carbon-D-Glucose)" %(step)
        solutionList = self._findGapfillSolution(stdMinimalModel)
        if self._isObjectMissing('Model', stdMinimalIntModel):
            print '   Integrating Carbon-D-Glucose media solutions...'
            self._integrateSolutions(stdMinimalModel, stdMinimalIntModel, solutionList, None)
        else:
            print '   Found Integrated Carbon-D-Glucose media solution %s/%s' %(self.args.workspace, stdMinimalIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdMinimalFba = "%s.model.std.int.minimal.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth on Carbon-D-Glucose" %(step)
        if self._isObjectMissing('FBA', stdMinimalFba):
            print '   Running FBA on the integrated model...'
            self._runFBA(stdMinimalIntModel, stdMinimalFba, media="Carbon-D-Glucose")
        else:
            print '   Found existing FBA solution %s/%s' %(self.args.workspace, stdMinimalFba)
        obj = self._getObjectiveValue(stdMinimalFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Network-based gapfilled model did not grow on minimal media after gapfilling'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        print '=== Completed Network-based Gap Fill Workflow ==='

        return

    def runProbabilistic(self):
        print '=== Started Likelihood-based Gap Fill Workflow ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsMeta = self._checkIfWorkspaceExists(self.args.workspace)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s' %(self.args.workspace, self.args.genome)
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
            rxnprobsMeta = self._runCalculate(probanno, rxnprobs)
        else:
           print '  Found reaction probabilities %s/%s' %(self.args.workspace, rxnprobs)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probModel = '%s.model.pa' %(self.args.genome)
        print '+++ Step %d: Run likelihood-based gap fill on complete media' %(step)
        if self._isObjectMissing('Model', probModel):
            print '  Submitting job and saving likelihood-based gap fill model to %s/%s ...' %(self.args.workspace, probModel)
            self._gapfill(draftModel, probModel, rxnprobs)
        else:
            print '  Found likelihood-based gap fill model %s/%s' %(self.args.workspace, probModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find unintegrated likelihood-based gap fill solutions' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, probModel)
        solutionList = self._findGapfillSolution(probModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIntModel = "%s.model.pa.int" %(self.args.genome)
        print "+++ Step %d: Integrate likelihood-based gap fill solution 1 on complete media (NOTE - you should check that the solution is optimal)" %(step)
        print '  Integrating solution %s into model %s/%s ...' %(solutionList[0], self.args.workspace, probIntModel)
        if self._isObjectMissing('Model', probIntModel):
            print '  Integrating probabilisitc gap fill solution 0 into model %s/%s ...' %(self.args.workspace, probIntModel)
            self._integrateSolutions(probModel, probIntModel, solutionList, rxnprobs)
        else:
            print '  Found integrated likelihood-based gap fill model %s/%s' %(self.args.workspace, probIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probCompleteFba = "%s.model.pa.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of likelihood-based gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('FBA', probCompleteFba):
            print '  Running fba and saving complete media likelihood-based gap fill FBA object to %s/%s' %(self.args.workspace, probCompleteFba)
            self._runFBA(probIntModel, probCompleteFba)
        else:
            print '  Found complete media likelihood-based gap fill FBA object %s/%s'  %(self.args.workspace, probCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        obj = self._getObjectiveValue(probCompleteFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Likelihood-based gapfilled model did not grow on complete media'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        step += 1
        probMinimalModel = "%s.model.pa.int.minimal" %(self.args.genome)
        print "+++ Step %d: Gapfill to minimal media (Carbon-D-Glucose)" %(step)
        if self._isObjectMissing('Model', probMinimalModel):
            print '   Gapfilling to Minimal media (Carbon-D-Glucose)'
            self._gapfill(probIntModel, probMinimalModel, rxnprobs, media='Carbon-D-Glucose')
        else:
            print '   Found carbon-D-glucose gapfilled model %s/%s' %(self.args.workspace, probMinimalModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probMinimalIntModel = "%s.model.pa.int.minimal.int" %(self.args.genome)
        print "+++ Step %d: Integrate gapfill solution to minimal media (Carbon-D-Glucose)" %(step)
        solutionList = self._findGapfillSolution(probMinimalModel)
        if self._isObjectMissing('Model', probMinimalIntModel):
            print '   Integrating Carbon-D-Glucose media solutions...'
            self._integrateSolutions(probMinimalModel, probMinimalIntModel, solutionList, rxnprobs)
        else:
            print '   Found Integrated Carbon-D-Glucose media solution %s/%s' %(self.args.workspace, probMinimalIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probMinimalFba = "%s.model.pa.int.minimal.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth on Carbon-D-Glucose" %(step)
        if self._isObjectMissing('FBA', probMinimalFba):
            print '   Running FBA on the integrated model...'
            self._runFBA(probMinimalIntModel, probMinimalFba, media="Carbon-D-Glucose")
        else:
            print '   Found existing FBA solution %s/%s' %(self.args.workspace, probMinimalFba)
        obj = self._getObjectiveValue(probMinimalFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Likelihood-based gapfilled model did not grow on minimal media after gapfilling'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        print '=== Completed Likelihood-based Gap Fill Workflow ==='

        return

    def runIterative(self):
        print '=== Started Network-based Iterative Gap Fill Workflow ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsmeta = self._checkIfWorkspaceExists(self.args.workspace)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s' %(self.args.workspace, self.args.genome)
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

        # NOTE - This is deliberately the same as for non-iterative gapfill.
        # We start with the integrated model from regular gapfill when doing iterative gapfill
        # because it is way faster and probably also more accurate to work around the MFA toolkit's
        # integration bug...
        step += 1
        stdModel = '%s.model.std' %(self.args.genome)
        print '+++ Step %d: Run network-based gap fill on complete media' %(step)
        if self._isObjectMissing('Model', stdModel):
            print '  Submitting job and saving network-based gap fill model to %s/%s ...' %(self.args.workspace, stdModel)
            self._gapfill(draftModel, stdModel, None)
        else:
            print '  Found network-based gap fill model %s/%s' %(self.args.workspace, stdModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find network-based gap fill unintegrated solutions (we only will integrate the first solution)' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, stdModel)
        solutionList = self._findGapfillSolution(stdModel)
        print '  Found unintegrated solution %s' %(solutionList[0])
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIntModel = "%s.model.std.int" %(self.args.genome)
        print '+++ Step %d: Integrate network-based gap fill solution 1 on complete media (NOTE - you should check that the solution is optimal)' %(step)
        if self._isObjectMissing('Model', stdIntModel):
            print '  Integrating solution %s into model %s/%s ...' %(solutionList[0], self.args.workspace, stdIntModel)
            self._integrateSolutions(stdModel, stdIntModel, solutionList, None)
        else:
            print '  Found integrated network-based gap fill model %s/%s' %(self.args.workspace, stdIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        
        # Note - iterative gapfilling automatically integrates solutions now.
        step += 1
        stdIterativeModel = '%s.model.std.iterative.int' %(self.args.genome)
        print '+++ Step %d: Run network-based iterative gap fill on complete media [note - iterative gapfill automatically will integrate the solutions]' %(step)
        if self._isObjectMissing('Model', stdIterativeModel):
            print '  Submitting job and saving network-based iterative gap fill model to %s/%s ...' %(self.args.workspace, stdIterativeModel)
            self._gapfill(stdIntModel, stdIterativeModel, None, iterative=True)
        else:
            print '  Found network-based iterative gap fill model %s/%s' %(self.args.workspace, stdIterativeModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        stdIterativeIntModel = stdIterativeModel

        step += 1
        stdIterativeIntSensitivity = '%s.model.std.iterative.int.sensitivity' %(self.args.genome)
        print '+++ Step %d: Check sensitivity of removing gap filled reactions...' %(step)
        if self._isObjectMissing('RxnSensitivity', stdIterativeIntSensitivity):
            print '  Submitting reaction sensitivity job and saving results to %s/%s ...' %(self.args.workspace, stdIterativeIntSensitivity)
            self._runReactionSensitivity(stdIterativeIntModel, stdIterativeIntSensitivity, rxnprobs=None)
        else:
            print '  Found rxn sensntivity object %s/%s' %(self.args.workspace, stdIterativeIntSensitivity)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIterativeIntModelFiltered = '%s.model.std.iterative.int.filtered' %(self.args.genome)
        print '+++ Step %d: Remove gapfilled reactions that are not necessary according to reaction sensitivity analysis' %(step)
        if self._isObjectMissing('Model', stdIterativeIntModelFiltered):
            print ' Deleteting unnecessary reactions from the model according to rxn sensitivity analysis... '
            self._runReactionDeletion(stdIterativeIntModel, stdIterativeIntModelFiltered, stdIterativeIntSensitivity)
        else:
            print '  Found filtered model %s/%s' %(self.args.workspace, stdIterativeIntModelFiltered)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))


        step += 1
        stdIterativeCompleteFba = "%s.model.std.iterative.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of network-based iterative gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('FBA', stdIterativeCompleteFba):
            print '  Running fba and saving complete media network-based iterative gap fill FBA object to %s/%s' %(self.args.workspace, stdIterativeCompleteFba)
            self._runFBA(stdIterativeIntModelFiltered, stdIterativeCompleteFba)
        else:
            print '  Found complete media network-based iterative gap fill FBA object %s/%s'  %(self.args.workspace, stdIterativeCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        obj = self._getObjectiveValue(stdIterativeCompleteFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Network-based iterative gapfilled model did not grow on complete media'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        step += 1
        stdIterativeMinimalModel = "%s.model.std.iterative.int.minimal" %(self.args.genome)
        print "+++ Step %d: Gapfill to minimal media (Carbon-D-Glucose)" %(step)
        if self._isObjectMissing('Model', stdIterativeMinimalModel):
            print '   Gapfilling to Minimal media (Carbon-D-Glucose)'
            self._gapfill(stdIterativeIntModelFiltered, stdIterativeMinimalModel, None, media='Carbon-D-Glucose')
        else:
            print '   Found carbon-D-glucose gapfilled model %s/%s' %(self.args.workspace, stdIterativeMinimalModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIterativeMinimalIntModel = "%s.model.std.iterative.int.minimal.int" %(self.args.genome)
        print "+++ Step %d: Integrate gapfill solution to minimal media (Carbon-D-Glucose)" %(step)
        solutionList = self._findGapfillSolution(stdIterativeMinimalModel)
        if self._isObjectMissing('Model', stdIterativeMinimalIntModel):
            print '   Integrating Carbon-D-Glucose media solutions...'
            self._integrateSolutions(stdIterativeMinimalModel, stdIterativeMinimalIntModel, solutionList, None)
        else:
            print '   Found Integrated Carbon-D-Glucose media solution %s/%s' %(self.args.workspace, stdIterativeMinimalIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        stdIterativeMinimalFba = "%s.model.std.iterative.int.minimal.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth on Carbon-D-Glucose" %(step)
        if self._isObjectMissing('FBA', stdIterativeMinimalFba):
            print '   Running FBA on the integrated model...'
            self._runFBA(stdIterativeMinimalIntModel, stdIterativeMinimalFba, media="Carbon-D-Glucose")
        else:
            print '   Found existing FBA solution %s/%s' %(self.args.workspace, stdIterativeMinimalFba)
        obj = self._getObjectiveValue(stdIterativeMinimalFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Network-based iterative gapfilled model did not grow on minimal media after gap filling.'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        print '=== Completed Network-based Iterative Gap Fill Workflow ==='

        return

    def runIterativeProbabilistic(self):
        print '=== Iterative Gap Fill Workflow (with likelihoods) ==='

        step = 0
        print '+++ Step %d: Initialize' %(step)
        print '  Probabilistic annotation service url is %s' %(self.args.paurl)
        print '  Workspace service url is %s' %(self.args.wsurl)
        print '  FBA model service url is %s' %(self.args.fbaurl)
        print '  Checking workspace %s ...' %(self.args.workspace)
        wsmeta = self._checkIfWorkspaceExists(self.args.workspace)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Load genome from the specified source' %(step)
        if self._isObjectMissing('Genome', self.args.genome):
            print '  Loading genome to %s/%s ...' %(self.args.workspace, self.args.genome)
            genomeMeta = self._loadGenome()
        else:
            print '  Found genome %s/%s' %(self.args.workspace, self.args.genome)
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
            rxnprobsMeta = self._runCalculate(probanno, rxnprobs)
        else:
           print '  Found reaction probabilities %s/%s' %(self.args.workspace, rxnprobs)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        # This is deliberately the same as for non-iterative gapfill. We want to start with the integrated solution
        # to biomass to work aroudn MFA toolkit's integration bug. It makes iterative gapfill SO much faster.
        step += 1
        probModel = '%s.model.pa' %(self.args.genome)
        print '+++ Step %d: Run likelihood-based gap fill on complete media' %(step)
        if self._isObjectMissing('Model', probModel):
            print '  Submitting job and saving likelihood-based gap fill model to %s/%s ...' %(self.args.workspace, probModel)
            self._gapfill(draftModel, probModel, rxnprobs)
        else:
            print '  Found likelihood-based gap fill model %s/%s' %(self.args.workspace, probModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Find unintegrated solutions for likelihood-based gap fill' %(step)
        print '  Checking gap fill model %s/%s ...' %(self.args.workspace, probModel)
        solutionList = self._findGapfillSolution(probModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIntModel = "%s.model.pa.int" %(self.args.genome)
        print "+++ Step %d: Integrate likelihood-based gap fill solution 1 on complete media (NOTE - you should check that the solution is optimal)" %(step)
        print '  Integrating solution %s into model %s/%s ...' %(solutionList[0], self.args.workspace, probIntModel)
        if self._isObjectMissing('Model', probIntModel):
            print '  Integrating likelihood-based gap fill solution 1 into model %s/%s ...' %(self.args.workspace, probIntModel)
            self._integrateSolutions(probModel, probIntModel, solutionList, rxnprobs)
        else:
            print '  Found integrated likelihood-based gap fill model %s/%s' %(self.args.workspace, probIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        #### Above this line, everything about this is exactly the same as for probanno without iterative gapfill

        step += 1
        probIterativeModel = '%s.model.pa.iterative.int' %(self.args.genome)
        print '+++ Step %d: Run likelihood-based iterative gap fill on complete media [ note - iterative gapfill will automatically integrate the solutions]' %(step)
        if self._isObjectMissing('Model', probIterativeModel):
            print '  Submitting job and saving likelihood-based iterative gap fill model to %s/%s ...' %(self.args.workspace, probIterativeModel)
            self._gapfill(probIntModel, probIterativeModel, rxnprobs, iterative=True)
        else:
            print '  Found likelihood-based gap fill model %s/%s' %(self.args.workspace, probIterativeModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        probIterativeIntModel = probIterativeModel

        step += 1
        probIterativeIntSensitivity = '%s.model.pa.iterative.int.sensitivity' %(self.args.genome)
        print '+++ Step %d: Order gapfilled reactions by probability and iteratively check the sensitivity of removing them' %(step)
        if self._isObjectMissing('RxnSensitivity', probIterativeIntSensitivity):
            print '  Submitting reaction sensitivity job and saving results to %s/%s ...' %(self.args.workspace, probIterativeIntSensitivity)
            self._runReactionSensitivity(probIterativeIntModel, probIterativeIntSensitivity, rxnprobs=rxnprobs)
        else:
            print '  Found rxn sensntivity object %s/%s' %(self.args.workspace, probIterativeIntSensitivity)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIterativeIntModelFiltered = '%s.model.pa.iterative.int.filtered' %(self.args.genome)
        print '+++ Step %d: Remove gapfilled reactions that are not necessary according to reaction sensitivity analysis' %(step)
        if self._isObjectMissing('Model', probIterativeIntModelFiltered):
            print ' Deleteting unnecessary reactions from the model according to rxn sensitivity analysis... '
            self._runReactionDeletion( probIterativeIntModel, probIterativeIntModelFiltered, probIterativeIntSensitivity)
        else:
            print '  Found filtered model %s/%s' %(self.args.workspace, probIterativeIntModelFiltered)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIterativeCompleteFba = "%s.model.pa.iterative.int.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth of likelihood-based iterative gap fill model on complete media (all available transporters to the cell are turned on)" %(step)
        if self._isObjectMissing('FBA', probIterativeCompleteFba):
            print '  Running fba and saving complete media likelihood-based iterative gap fill FBA object to %s/%s' %(self.args.workspace, probIterativeCompleteFba)
            self._runFBA(probIterativeIntModelFiltered, probIterativeCompleteFba)
        else:
            print '  Found complete media likelihood-based iterative gap fill FBA object %s/%s'  %(self.args.workspace, probIterativeCompleteFba)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        obj = self._getObjectiveValue(probIterativeCompleteFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Likelihood-based iterative gapfilled model did not grow on complete media'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        step += 1
        probIterativeMinimalModel = "%s.model.pa.iterative.int.minimal" %(self.args.genome)
        print "+++ Step %d: Gapfill to minimal media (Carbon-D-Glucose)" %(step)
        if self._isObjectMissing('Model', probIterativeMinimalModel):
            print '   Gapfilling to Minimal media (Carbon-D-Glucose)'
            self._gapfill(probIterativeIntModelFiltered, probIterativeMinimalModel, rxnprobs, media='Carbon-D-Glucose')
        else:
            print '   Found carbon-D-glucose gapfilled model %s/%s' %(self.args.workspace, probIterativeMinimalModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIterativeMinimalIntModel = "%s.model.pa.iterative.int.minimal.int" %(self.args.genome)
        print "+++ Step %d: Integrate gapfill solution to minimal media (Carbon-D-Glucose)" %(step)
        solutionList = self._findGapfillSolution(probIterativeMinimalModel)
        if self._isObjectMissing('Model', probIterativeMinimalIntModel):
            print '   Integrating Carbon-D-Glucose media solutions...'
            self._integrateSolutions(probIterativeMinimalModel, probIterativeMinimalIntModel, solutionList, rxnprobs)
        else:
            print '   Found Integrated Carbon-D-Glucose media solution %s/%s' %(self.args.workspace, probIterativeMinimalIntModel)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        probIterativeMinimalFba = "%s.model.pa.iterative.int.minimal.fba" %(self.args.genome)
        print "+++ Step %d: Check for growth on Carbon-D-Glucose" %(step)
        if self._isObjectMissing('FBA', probIterativeMinimalFba):
            print '   Running FBA on the integrated model...'
            self._runFBA(probIterativeMinimalIntModel, probIterativeMinimalFba, media="Carbon-D-Glucose")
        else:
            print '   Found existing FBA solution %s/%s' %(self.args.workspace, probIterativeMinimalFba)
        obj = self._getObjectiveValue(probIterativeMinimalFba)
        if obj is None or float(obj) < 1E-5:
            print '   [ERROR] Likelihood-based iterative gapfilled model did not grow on minimal media after gapfilling'
            raise NoGrowthError('Growth rate was 0 or unable to calculate')

        print '=== Completed Likelihood-based Iterative Gap Fill Workflow ==='
        
        return

    def runKnockouts(self):
        print "=== Started Knockout Workflow ==="

        # We need to turn this setting off in a few places
        oldforce = self.args.force

        step = 0
        print '+++ Step %d: Search for integrated gap fill models to integrate knockout data +++' %(step)
        integratedModelNames = [ "%s.model.std.int.minimal.int" %(self.args.genome),
                                 "%s.model.pa.int.minimal.int" %(self.args.genome),
                                 "%s.model.std.iterative.int.minimal.int" %(self.args.genome),
                                 "%s.model.pa.iterative.int.minimal.int" %(self.args.genome) ]

        self.args.force = False
        foundModels = []
        for model in integratedModelNames:
            if self._isObjectMissing("Model", model):
                continue
            else:
                foundModels.append(model)
        if len(foundModels) == 0:
            print "  [ERROR] No gap fill models found for which to integrate knockout data. Run the workflow with a gapfill flag first."
            raise NoDataError("No integrated models found")
        else:
            for model in foundModels:
                print '  Found integrated gap fill model %s/%s' %(self.args.workspace, model)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        self.args.force = oldforce

        step += 1
        print '+++ Step %d: Find all media conditions in phenotype set +++' %(step)
        print '  Checking phenotype set %s/%s...' %(self.args.workspace, self.args.knockout)
        mediaTuples = self._getMediaIdsFromPhenotypeSet(self.args.knockout, workspace=self.args.knockoutws)
        print '  Found %d media conditions in phenotype set' %(len(mediaTuples))
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        step += 1
        print '+++ Step %d: Run gap fill and integrate solutions for each media found in the phenotype set +++' %(step)
        for model in foundModels:
            for media in mediaTuples:
                mediaGapfilledModel = model + ".%s" %(media[0])
                mediaGapfilledIntModel = mediaGapfilledModel + ".int"

                # Does the integrated model ALREADY GROW on this media?
                substep = 1
                print '+++ Step %d.%d: Check for growth of %s model on %s media +++' %(step, substep, model, media[0])
                mediaTestFba = model + ".%s.testgrowth.fba" %(media[0])
                if self._isObjectMissing("FBA", mediaTestFba):
                    self._runFBA(model, mediaTestFba, media=media[0], mediaws=media[1])
                else:
                    print "  Found FBA solution object %s/%s" %(self.args.workspace, mediaTestFba)
                obj = self._getObjectiveValue(mediaTestFba)
                print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

                if obj is None or float(obj) < 1E-5:
                    substep += 1
                    print '+++ Step %d.%d: Determine type of gap fill for media %s +++' %(step, substep, media[0])
                    # First we need to see if the gapfill was a standard (std) or probanno (pa) gapfill.
                    self.args.force = False
                    isProbanno = False
                    rxnprobs = None
                    if ".pa." in model:
                        isProbanno = True
                        print '  Type of gap fill is likelihood-based'
                        rxnprobs = '%s.rxnprobs' %(self.args.genome)
                        print '  Looking for RxnProbs object %s/%s...' %(self.args.workspace, rxnprobs)
                        if self._isObjectMissing('RxnProbs', rxnprobs):
                            print '  [ERROR] RxnProbs object was not found. Try running one of the likelihood-based gapfilling options of this workflow again first. '
                            raise NoErrorData("No rxnprobs object found")
                        print '  Found RxnProbs object'
                    else:
                        print '  Type of gap fill is network-based'
                    self.args.force = oldforce
                    print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
                     
                    substep += 1
                    print '+++ Step %d.%d: Run gap fill on %s media +++' %(step, substep, media[0])
                    if self._isObjectMissing('Model', mediaGapfilledModel):
                        print '  Submitting job and saving gap fill model to %s/%s...' %(self.args.workspace, mediaGapfilledModel)
                        self._gapfill(model, mediaGapfilledModel, rxnprobs, iterative=False, media=media[0], mediaws=media[1])
                    else:
                        print '  Found gap fill model %s/%s' %(self.args.workspace, mediaGapfilledModel)
                    print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))             
    
                    substep += 1
                    print '+++ Step %d.%d: Find and integrate solution on %s media +++' %(step, substep, media[0])
                    solutionList = self._findGapfillSolution(mediaGapfilledModel)
                    if self._isObjectMissing('Model', mediaGapfilledIntModel):
                        print '  Integrating solution %s into model %s/%s' %(solutionList[0], self.args.workspace, mediaGapfilledIntModel)
                        self._integrateSolutions(mediaGapfilledModel, mediaGapfilledIntModel, solutionList, rxnprobs)
                    else:
                        print '  Found integrated model %s/%s' %(self.args.workspace, mediaGapfilledIntModel)
                    print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
    
                    substep += 1
                    print '+++ Step %d.%d: Check for growth in integrated model on media %s +++' %(step, substep, media[0])
                    finalFbaObject = "%s.fba" %(mediaGapfilledIntModel)
                    if self._isObjectMissing('FBA', finalFbaObject):
                        print '  Running fba and saving FBA object to %s/%s' %(self.args.workspace, finalFbaObject)
                        self._runFBA(mediaGapfilledIntModel, finalFbaObject, media=media[0], mediaws=media[1])
                    else:
                        print '  Found FBA object %s/%s' %(self.args.workspace, finalFbaObject)
                    print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
    
                    obj = self._getObjectiveValue(finalFbaObject)
                    if obj is None or obj < 1E-5:
                        raise NoGrowthError("Did not get growth after gapfilling to knockout media %s" %(media[0]))
                else:
                    substep += 1
                    print "+++ Step %d.%d: No gapfilling needed - model already grows on %s media. +++" %(step, substep, media[0])
                    mediaGapfilledModel = model
                    mediaGapfilledIntModel = model
                    print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

            # Do the simulation - get simulation set (note - for knockouts we don't want to necessarily add transporters,
            # they should already be there from the gapfill we did above)
            substep += 1
            print '+++ Step %d.%d: simulate phenotype +++' %(step, substep)
            knockoutSimulation = "%s.knockoutsim" %(model)
            if self._isObjectMissing('PhenotypeSimulationSet', knockoutSimulation):
                print '  Running phenotype simulation and saving to %s/%s' %(self.args.workspace, knockoutSimulation)
                self._simulatePhenotype(mediaGapfilledIntModel, self.args.knockout, self.args.knockoutws, knockoutSimulation, positive_transporters = 0, all_transporters = 0)
            else:
                print '  Found phenotype simulation set %s/%s' %(self.args.workspace, knockoutSimulation)
            print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))

        print '=== Completed Knockout Workflow ==='
        
        return
                
    def runBiolog(self):

        print "=== Running Biolog workflow ==="

        # We need to turn this setting off in a few places
        oldforce = self.args.force

        if self.args.positiveTransportersOnly:
            positive_transporters = 1
            all_transporters = 0
        else:
            positive_transporters = 0
            all_transporters = 1

        step = 0
        print '+++ Step %d: Search for integrated gap fill models (to Carbon-D-Glucose) to use to test biolog data +++' %(step)
        integratedModelNames = [ "%s.model.std.int.minimal.int" %(self.args.genome),
                                 "%s.model.pa.int.minimal.int" %(self.args.genome),
                                 "%s.model.std.iterative.int.minimal.int" %(self.args.genome),
                                 "%s.model.pa.iterative.int.minimal.int" %(self.args.genome) ]

        self.args.force = False
        foundModels = []
        for model in integratedModelNames:
            if self._isObjectMissing("Model", model):
                continue
            else:
                foundModels.append(model)
        if len(foundModels) == 0:
            print "  [ERROR] No gap fill models found for which to integrate knockout data. Run the workflow with a gapfill flag first."
            raise NoDataError("No integrated models found")
        else:
            for model in foundModels:
                print '  Found integrated gap fill model %s/%s' %(self.args.workspace, model)
        print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
        self.args.force = oldforce

        step += 1
        substep = 1
        print '+++ Step %d: Run biolog data simulation on each of the found models. +++' %(step)
        for model in foundModels:
            print "+++ Step %d/%d: Simulate biolog data on model %s/%s" %(step, substep, self.args.workspace, model)
            biologSimulation = "%s.biologsim" %(model)
            if self._isObjectMissing("PhenotypeSimulationSet", biologSimulation):
                print "   Running simulation..."
                self._simulatePhenotype(model, self.args.biolog, self.args.biologws, biologSimulation, positive_transporters = positive_transporters, all_transporters = all_transporters)
            else:
                print "   Found existing phenotype simulation object %s" %(biologSimulation)
            print '  [OK] %s' %(time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime()))
            substep += 1

        print "=== Completed biolog simulation workflow ==="

        return

if __name__ == "__main__":
    # Parse options.
    description=''' Workflow.py

    DESCRIPTION:
           Performs all of the steps in the workflows described in the manuscript. Each of the following four flags indicates a different workflow:

           --network :             Perform network-based targeted gap filling (target is biomass)
           --likelihood:           Perform likelihood-based targeted gap filling (target is biomass)
           --iterative:            Perform iterative network-based gap filling workflow
           --iterative-likelihood: Perform iterative likelihood-based gap filling workflow

           The only required inputs are genome and --workspace. The genome is automatically loaded from the appropriate database (by defualt,
           this database is the PubSEED: http://pubseed.theseed.org) and then the chosen workflow (above) is run.

           Phenotype data consistency is tested separately because it has to be done with extra external data. 

           If you have KNOCKOUT data: Put it in a PhenotypeSet object and run this script with:
               --knockout (name of PhenotypeSet object)
               --knockoutws (Optional: Name of workspace in which PhenotypeSet object is found)

           The knockout data workflow differs from Biolog in that we will make sure the model is capable of growing on media on which the knockouts
           were performed before testing the effects of knockouts. If it is not we run a gap filling to fix it.

           If you have BIOLOG data: Put it in a PhenotypeSet object and run 
               --biologdata (name of PhenotypeSet object)
               --biologdataws (Optional: Name of workspace in which PhenotypeSet object is found)

           By defualt when simulating Biolog data we try to add transporters for every media condition present in the 
           biolog dataset. If you want to try to add transporters for ONLY positive-growth conditions (a form of data fitting),
           use --positiveTransportersOnly.

           In either case you have to have run at least one of the other workflows first (or concurrently) to generate models before running the
           phenotype simulations.

           By default the URLs point at production servers. If you are a KBase developer and want to test effects of using a modified service, you
           can change the URLs for these services.

    '''

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='Workflow', description=description)
    parser.add_argument('genome', help='ID of genome that matches the source (e.g. 83333.1 for SEED, kb|g.0 for kbase)', action='store', default=None)
    parser.add_argument('--genome_source', help='Source for genome. Valid values include "kbase", "seed", "rast". (Default is seed)', action='store', dest='source', default='seed')
    parser.add_argument('-w', '--workspace', help='workspace for storing objects', action='store', dest='workspace')
    parser.add_argument('--force', help='Force rebuilding of all objects', action='store_true', dest='force', default=False)
    parser.add_argument('--network', help='Run network-based gap fill workflow', action='store_true', dest='standard', default=False)
    parser.add_argument('--likelihood', help='Run likelihood-based gap fill workflow', action='store_true', dest='prob', default=False)
    parser.add_argument('--iterative', help='Run network-based iterative gap fill workflow ', action='store_true', dest='iterative', default=False)
    parser.add_argument('--iterative-likelihood', help='Run likelihood-based iterative gap fill workflow ', action='store_true', dest='iterativeprob', default=False)
    parser.add_argument('--num-solutions', help='Number of solutions to find for gap fill (ignored for iterative gap filling)', action='store', dest='numsolutions', type=int, default=3)
    parser.add_argument('--ws-url', help='URL for workspace service', action='store', dest='wsurl', default='http://kbase.us/services/ws')
    parser.add_argument('--job-url', help='URL for old workspace service (needed for job support)', action='store', dest='joburl', default='http://kbase.us/services/workspace')
    parser.add_argument('--fba-url', help='URL for fba model service', action='store', dest='fbaurl', default='http://kbase.us/services/fba_model_services')
    parser.add_argument('--pa-url', help='URL for probabilistic annotation service', action='store', dest='paurl', default='http://kbase.us/services/probabilistic_annotation')
    parser.add_argument('--knockout', help='OPTIONAL. Provide a knockout data PhenotypeSet and we will create a new model gapfilled to each media in it.', action='store', dest='knockout', default=None)
    parser.add_argument('--knockoutws', help='OPTIONAL. Workspace for provided knockout data PhenotypeSet (default is the same as the current workspace)', action='store', dest='knockoutws', default=None)
    parser.add_argument('--biologdata', help='OPTIONAL. Biolog PhenotypeSet object', action='store', dest='biolog', default=None)
    parser.add_argument('--biologdataws', help='OPTIONAL. Workspace for provided biolog data PhenotypeSet (default is the same as the current workspace)', action='store', dest='biologws', default=None)
    parser.add_argument('--positiveTransportersOnly', help='OPTIONAL. Add transporters for only media that produced growth in a PhenotypeSet before simulating. Only relevant for biolog data. Default is to add transporters for ALL media.', action='store_true', default=False)
    parser.add_argument('--maxtime', 
                        help='OPTIONAL. Maximum amount of time to wait for a job to finish (by default the maximum is 2 hours for targeted gapfill jobs and 4 days for iterative gapfill)',
                        action='store', default=None)
    parser.add_argument('--ujsurl', help='URL for user and job state service (needed to check probanno jobs)', action='store', dest='ujsurl', default='https://kbase.us/services/userandjobstate/')
    args = parser.parse_args()

    if args.knockoutws is None:
        args.knockoutws = args.workspace

    if args.biologws is None:
        args.biologws = args.workspace
    
    try:
        if args.force:
            print '[WARNING] All objects will be rebuilt because --force option was specified.'
        workflow = Workflow(args)
        if workflow.args.standard:
            oldmax = workflow.args.maxtime
            if workflow.args.maxtime is None:
                workflow.args.maxtime = 86400
            workflow.runStandard()
            workflow.args.maxtime = oldmax
        if workflow.args.prob:
            oldmax = workflow.args.maxtime
            if workflow.args.maxtime is None:
                workflow.args.maxtime = 86400
            workflow.runProbabilistic()
            workflow.args.maxtime = oldmax
        if workflow.args.iterative:
            oldmax = workflow.args.maxtime
            if workflow.args.maxtime is None:
                # 4 days
                workflow.args.maxtime = 345600
            workflow.runIterative()
            workflow.args.maxtime = oldmax
        if workflow.args.iterativeprob:
            oldmax = workflow.args.maxtime
            if workflow.args.maxtime is None:
                workflow.args.maxtime = 345600
            workflow.runIterativeProbabilistic()
            workflow.args.maxtime = oldmax
        if workflow.args.knockout is not None:
            if workflow.args.maxtime is None:
                workflow.args.maxtime = 86400
            workflow.runKnockouts()
        if workflow.args.biolog is not None:
            workflow.runBiolog()

    except Exception as e:
        print '  [ERROR]'
        traceback.print_exc(file=sys.stderr)
        exit(1)

    exit(0)


        
