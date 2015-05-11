
import subprocess
import sys
import os
import shutil
import traceback
import time
import math
import re
import tempfile
from biokbase.probabilistic_annotation.Helpers import make_object_identity, make_job_directory, timestamp, ProbAnnoType, RxnProbsType, WrongVersionError, ServiceVersion, ServiceName, get_config
from biokbase.probabilistic_annotation.DataParser import DataParser
from biokbase.userandjobstate.client import UserAndJobState
from biokbase.workspace.client import Workspace
from biokbase.cdmi.client import CDMI_EntityAPI
from biokbase import log
from urllib2 import HTTPError

# Exception thrown when no features are found in Genome object
class NoFeaturesError(Exception):
    pass

# Exception thrown when blast command failed
class BlastError(Exception):
    pass

# Exception thrown when there is an invalid number calculating likelihoods
class BadLikelihoodError(Exception):
    pass

# Exception thrown when a target id is not found in rolestring dictionary
class NoTargetIdError(Exception):
    pass

# Exception thrown when there are no gene IDs in Genome object
class NoGeneIdsError(Exception):
    pass

# Exception thrown when role not found in roleToTotalProb dictionary
class RoleNotFoundEror(Exception):
    pass

''' Worker that implements probabilistic annotation algorithm. '''

class ProbabilisticAnnotationWorker:

    def runAnnotate(self, jobID, input):
        ''' Calculate annotation likelihoods and store them in a ProbAnno typed object.

            A ProbAnno typed object is created in four steps: (1) extract amino acid
            sequences from a Genome typed object to a fasta file, (2) run a protein
            search using the amino acid sequences against the protein database,
            (3) calculate annotation likelihood scores for each roleset implied by the
            functions of proteins in subsystems, and (4) save the likelihood scores
            to a ProbAnno typed object.

            @param jobID: Job ID from user and job state service
            @param input: Input parameters for annotate() function
            @return Nothing (although job is marked as complete)
        '''

        # Assuming for now that a logger can be setup.
        if input['verbose'] and self.logger.get_log_level() < log.DEBUG:
            self.logger.set_log_level(log.DEBUG)

        # Get the static database files.  If the files do not exist and they are downloaded
        # from Shock, it can take a few minutes before they are ready.
        testDataPath = os.path.join(os.environ['KB_TOP'], 'services', 'testdata')
        self.dataParser.getDatabaseFiles(self.logger, testDataPath)

        # Initial job status.
        status = None

        try:
            # Make sure the database files are available.
            self.dataParser.checkIfDatabaseFilesExist()

            # Make sure the job directory exists.
            workFolder = make_job_directory(self.config['work_folder_path'], jobID)

            # Create a user and job state client and authenticate as the user.
            ujsClient = UserAndJobState(self.config['userandjobstate_url'], token=self.ctx['token'])
    
            # Get the Genome object from the specified workspace.
            try:
                ujsClient.update_job_progress(jobID, self.ctx['token'], 'getting genome object', 1, timestamp(3600))
            except:
                pass
            wsClient = Workspace(self.config['workspace_url'], token=self.ctx['token'])
            genomeObjectId = make_object_identity(input['genome_workspace'], input['genome'])
            objectList = wsClient.get_objects( [ genomeObjectId ] )
            genomeObject = objectList[0]
            
            # Convert Genome object to fasta file.
            try:
                ujsClient.update_job_progress(jobID, self.ctx['token'], 'converting Genome object to fasta file', 1, timestamp(3600))
            except:
                pass
            fastaFile = self._genomeToFasta(input, genomeObject, workFolder)
            
            # Run blast using the fasta file.
            try:
                ujsClient.update_job_progress(jobID, self.ctx['token'], 'running blast', 1, timestamp(3600))
            except:
                pass
            blastResultFile = self._runBlast(input, fastaFile, workFolder)
            
            # Calculate roleset probabilities.
            try:
                ujsClient.update_job_progress(jobID, self.ctx['token'], 'calculating roleset probabilities', 1, timestamp(300))
            except:
                pass
            rolestringTuples = self._rolesetProbabilitiesMarble(input, blastResultFile, workFolder)
            
            # Build ProbAnno object and store in the specified workspace.
            try:
                ujsClient.update_job_progress(jobID, self.ctx['token'], 'building ProbAnno object', 1, timestamp(120))
            except:
                pass
            output = self._buildProbAnnoObject(input, genomeObject, blastResultFile, rolestringTuples, workFolder, wsClient)

            # Mark the job as done.
            status = 'done'
            tb = None
            self._log(log.INFO, 'Job '+jobID+' finished for genome '+input['genome']+' to probanno '+input['probanno'])

        except:
            tb = traceback.format_exc()
            sys.stderr.write('[error] calculating annotation likelihoods failed\n'+tb)
            status = 'failed'
            self._log(log.ERR, 'Job '+jobID+' failed for genome '+input['genome']+' to probanno '+input['probanno'])
        
        # Mark the job as complete with the given status.
        ujsClient.complete_job(jobID, self.ctx['token'], status, tb, { })

        # Remove the temporary work directory.
        if self.logger.get_log_level() < log.DEBUG2 and status == 'done':
            try:
                shutil.rmtree(workFolder)
            except OSError:
                # For some reason deleting the directory was failing in production. Rather than have all jobs look like they failed
                # I catch and log the exception here (since the user still gets the same result if the directory remains intact)
                msg = 'Unable to delete temporary directory %s\n' %(workFolder)
                sys.stderr.write('[warning] '+msg)
                self._log(log.WARNING, msg)

        return
        
    def _genomeToFasta(self, input, genomeObject, workFolder):

        ''' Convert a Genome object into an amino-acid FASTA file (for BLAST purposes).

            @param input Dictionary of input parameters to annotate() function
            @param genomeObject Genome typed object from workspace
            @param workFolder Path to directory in which to store temporary files
            @return Path to fasta file with query proteins
            @raise NoFeaturesError when Genome object has no features
        '''

        # Make sure the Genome object has features.
        if 'features' not in genomeObject['data']:
            raise NoFeaturesError('The input Genome object %s/%s has no features. Did you forget to run annotate_genome?\n' %(input['genome_workspace'], input['genome']))
    
        # Run the list of features to build the fasta file.
        sys.stdout.write('Creating fasta file for genome %s\n' %(input['genome']))
        fastaFile = os.path.join(workFolder, '%s.faa' %(input['genome']))
        with open(fastaFile, 'w') as handle:
            numProteins = 0
            features = genomeObject['data']['features']
            for feature in features:
                # Not a protein-encoding gene
                if 'protein_translation' not in feature:
                    continue
                handle.write('>%s\n%s\n' %(feature['id'], feature['protein_translation']))
                numProteins += 1
        
        sys.stdout.write('Wrote %d protein sequences to "%s"\n' %(numProteins, fastaFile))
        self._log(log.DEBUG, 'Wrote %d protein sequences to "%s"' %(numProteins, fastaFile))
        return fastaFile
        
    def _runBlast(self, input, queryFile, workFolder):

        ''' A simplistic wrapper to search for the query proteins against the subsystem proteins.

            @param input Dictionary of input parameters to annotate() function
            @param queryFile Path to fasta file with query proteins
            @param workFolder Path to directory in which to store temporary files
            @return Path to output file from search program
            @raise BlastError when there is a problem running the search program
        '''

        # Generate path to output file.  Output format 6 is tab-delimited format.
        blastResultFile = os.path.join(workFolder, '%s.blastout' %(input['genome']))

        # Build the command based on the configured search program.
        if self.config['search_program'] == 'usearch':
            args = [ self.config['search_program_path'], '-ublast', queryFile,
                     '-db', self.dataParser.SearchFiles['protein_udb_file'],
                     '-evalue', self.config['search_program_evalue'],
                     '-accel', self.config['usearch_accel'],
                     '-threads', self.config['blast_threads'],
                     '-blast6out', blastResultFile ]
        else:
            args = [ self.config['search_program_path'], '-query', queryFile,
                     '-db', self.dataParser.DataFiles['protein_fasta_file'],
                     '-outfmt', '6', '-evalue', self.config['search_program_evalue'],
                     '-num_threads', self.config['blast_threads'],
                     '-out', blastResultFile ]

        # Run the command to search for proteins against subsystem proteins.
        cmd = ' '.join(args)
        message = 'Started protein search with command: '+cmd
        sys.stdout.write(message+'\n')
        self._log(log.INFO, message)
        try:
            proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (stdout, stderr) = proc.communicate()
            if proc.returncode < 0:
                message = '"%s" was terminated by signal %d' %(args[0], -proc.returncode)
                raise BlastError(message)
            else:
                if proc.returncode > 0:
                    details = '"%s" failed with return code %d\nCommand: "%s"\nStdout: "%s"\nStderr: "%s"' \
                        %(args[0], proc.returncode, cmd, stdout, stderr)
                    raise BlastError(details)
        except OSError as e:
            message = 'Failed to run "%s": %s' %(args[0], e.strerror)
            raise BlastError(message)
        sys.stdout.write('Finished protein search\n')

        return blastResultFile
    
    def _rolesetProbabilitiesMarble(self, input, blastResultFile, workFolder):

        ''' Calculate the probabilities of rolesets from the BLAST results.

            A roleset is each possible combination of roles implied by the functions
            of the proteins in subsystems.  The output is a dictionary keyed by
            query gene of lists of tuples where each tuple contains (1) roleset
            string, and (2) likelihood value.  The roleset string is a concatenation
            of all of the roles of a protein with a single function (order does
            not matter).
    
            @param input Dictionary of input parameters to annotate() function
            @param blastResultFile Path to output file from BLAST
            @param workFolder Path to directory in which to store temporary files
            @return Dictionary keyed by query gene of list of tuples with roleset and likelihood
            @raise BadLikelihoodError when there is math error calculating a likelihood
            @raise NoTargetIdError when target ID in search results is not found in rolestrings
        '''

        message = 'Started marble-picking on rolesets for genome %s' %(input['genome'])
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)
    
        # Read in the target roles (this function returns the roles as lists!)
        targetIdToRole, targetRoleToId = self.dataParser.readFidRoleFile(self.dataParser.DataFiles['otu_fid_role_file'])
    
        # Convert the lists of roles into "rolestrings" (sort the list so that order doesn't matter)
        # in order to deal with the case where some of the hits are multi-functional and others only have
        # a single function.
        targetIdToRoleString = dict()
        for target in targetIdToRole:
            stri = self.config['separator'].join(sorted(targetIdToRole[target]))
            targetIdToRoleString[target] = stri

        # Parse the output from BLAST which returns a dictionary keyed by query gene of a list
        # of tuples with target gene and score.
        # query --> [ (target1, score 1), (target 2, score 2), ... ]
        idToTargetList = self.dataParser.parseBlastOutput(blastResultFile)
    
        # This is a holder for all of our results which is a dictionary keyed by query gene
        # of a list of tuples with roleset and likelihood.
        # query -> [ (roleset1, likelihood_1), (roleset2, likelihood_2), ...]
        rolestringTuples = dict()

        # For each query gene we calculate the likelihood of each possible rolestring
        # See equation 2 in the paper ("Calculating annotation likelihoods" section).
        for query in idToTargetList:
            # First we need to know the maximum score for this gene.
            # I have no idea why but I'm pretty sure Python is silently turning the second
            # element of these tuples into strings.  That's why I turn them back to floats.
            maxscore = 0
            for tup in idToTargetList[query]:
                if float(tup[1]) > maxscore:
                    maxscore = float(tup[1])
    
            # Now we calculate the cumulative squared scores for each possible rolestring.
            # This along with pseudocount*maxscore is equivalent to multiplying all scores
            # by themselves and then dividing by the max score.
            # This is done to avoid some pathological cases and give more weight to higher-scoring hits
            # and not let much lower-scoring hits \ noise drown them out.
            # Build a dictionary keyed by rolestring of the sum of squares of the log-scores.
            rolestringToScore = dict()
            for tup in idToTargetList[query]:
                try:
                    rolestring = targetIdToRoleString[tup[0]]
                except KeyError:
                    message = 'Target id %s from search results file had no roles in rolestring dictionary' %(tup[0])
                    raise NoTargetIdError(message)
                if rolestring in rolestringToScore:
                    rolestringToScore[rolestring] += (float(tup[1]) ** 2)
                else:
                    rolestringToScore[rolestring] = (float(tup[1]) ** 2)
    
            # Calculate the likelihood that this gene has the given functional annotation.
            # Start with the denominator which is the sum of squares of the log-scores for
            # all possible rolestrings.
            denom = float(self.config['pseudo_count']) * maxscore
            for stri in rolestringToScore:
                denom += rolestringToScore[stri]
            if math.isnan(denom):
                message = 'Denominator in likelihood calculation for gene %s is NaN %f' %(query, denom)
                raise BadLikelihoodError(message)

            # The numerators are the sum of squares for each rolestring.
            # Calculate the likelihood for each rolestring and store in the output dictionary.
            for stri in rolestringToScore:
                p = rolestringToScore[stri] / denom
                if math.isnan(p):
                    message = 'Likelihood for rolestring %s in gene %s is NaN based on score %f' %(stri, query, rolestringToScore[stri])
                    raise BadLikelihoodError(message)
                if query in rolestringTuples:
                    rolestringTuples[query].append( (stri, p) )
                else:
                    rolestringTuples[query] = [ (stri, p) ]
    
        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            rolesetProbabilityFile = os.path.join(workFolder, '%s.rolesetprobs' %(input['genome']))
            with open(rolesetProbabilityFile, 'w') as handle:
                for query in rolestringTuples:
                    for tup in rolestringTuples[query]:
                        handle.write('%s\t%s\t%1.4f\n' %(query, tup[0], tup[1]))
            
        message = 'Finished marble-picking on rolesets for genome %s' %(input['genome'])
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)
        return rolestringTuples
            
    def _buildProbAnnoObject(self, input, genomeObject, blastResultFile, queryToRolesetProbs, workFolder, wsClient):

        ''' Create a ProbAnno typed object and save it to a workspace.

            The queryToRolesetProbs dictionary has this format: querygene -> [ (roleset, likelihood), ... ]
            The probabilistic annotation object adds fields for the probability of each role being linked to each gene.

            @param input Dictionary of input parameters to annotate() function
            @param genomeObject Genome typed object from workspace
            @param blastResultFile Path to output file from BLAST in tab-delimited format
            @param queryToRolesetProbs: Dictionary keyed by query protein of list of tuples with roleset and likelihood
            @param workFolder Path to directory in which to store temporary files
            @param wsClient Workspace client object
            @return metadata
            @raise NoGeneIdsError when there are no gene IDs in Genome object
        '''

        message = 'Started building ProbAnno object %s/%s for genome %s...' %(input['probanno_workspace'], input['probanno'], input['genome'])
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        # Read in the target roles (this function returns the roles as lists!)
        targetToRoles, rolesToTargets = self.dataParser.readFidRoleFile(self.dataParser.DataFiles['otu_fid_role_file'])
        targetToRoleSet = dict()
        for target in targetToRoles:
            stri = self.config['separator'].join(sorted(targetToRoles[target]))
            targetToRoleSet[target] = stri

        # This is a dictionary from query ID to (target, -log E-value) pairs.
        # We just use it to identify whether or not we actually hit anything in the db
        # when searching for the query gene.
        queryToTargetEvals = self.dataParser.parseBlastOutput(blastResultFile)
        
        # For each query ID:
        # 1. Identify their rolestring probabilities (these are the first and second elements of the tuple)
        # 2. Iterate over the target genes and identify those with each function (a list of these and their blast scores is
        #    the third element of the tuple) - should be able to set it up to only iterate once over the list.
        # 3. Add that tuple to the JSON file with the key "alternativeFunctions"
    
        # The Genome object data ["features"] is a list of dictionaries. We want to make our data structure and 
        # then add that to the dictionary.  I use the ii in range so I can edit the elements without changes being lost.
    
        objectData = dict()
        objectData['id'] = input['probanno']
        objectData['genome'] = input['genome']
        objectData['genome_workspace'] = input['genome_workspace'];
        objectData['roleset_probabilities'] = queryToRolesetProbs;
        objectData['skipped_features'] = []
        
        for ii in range(len(genomeObject['data']['features'])):
            feature = genomeObject['data']['features'][ii]
            if 'id' not in genomeObject['data']:
                raise NoGeneIdsError('No gene IDs found in input Genome object %s/%s (this should never happen)' %(input['genome_workspace'], input['genome']))
            queryid = feature['id']
    
            # This can happen if I couldn't find hits from that gene to anything in the database. In this case, I'll just skip it.
            # TODO Or should I make an empty object? I should ask Chris.
            if queryid not in queryToRolesetProbs or queryid not in queryToTargetEvals:
                objectData['skipped_features'].append(queryid)
                
        # Store the ProbAnno object in the specified workspace.
        objectMetaData = dict()
        objectMetaData['num_rolesets'] = len(objectData['roleset_probabilities'])
        objectMetaData['num_skipped_features'] = len(objectData['skipped_features'])
        objectProvData = dict()
        objectProvData['time'] = timestamp(0)
        objectProvData['service'] = os.environ.get('KB_SERVICE_NAME', ServiceName)
        objectProvData['service_ver'] = ServiceVersion
        objectProvData['method'] = 'annotate'
        objectProvData['method_params'] = input.items()
        objectProvData['input_ws_objects'] = [ '%s/%s/%d' %(genomeObject['info'][7], genomeObject['info'][1], genomeObject['info'][4]) ]
        objectSaveData = dict()
        objectSaveData['type'] = ProbAnnoType
        objectSaveData['name'] = input['probanno']
        objectSaveData['data'] = objectData
        objectSaveData['meta'] = objectMetaData
        objectSaveData['provenance'] = [ objectProvData ]
        retryCount = 3
        while retryCount > 0:
            try:
                objectInfo = wsClient.save_objects( { 'workspace': input['probanno_workspace'], 'objects': [ objectSaveData ] } )
                message = 'Saved ProbAnno object %s/%s' %(input['probanno_workspace'], input['probanno'])
                sys.stdout.write(message+'\n')
                self._log(log.DEBUG, message)
                return objectInfo[0]
            except HTTPError as e:
                # Hopefully this is just a temporary glitch, try again in a few seconds since we worked so hard to build the object.
                retryCount -= 1
                message = 'HTTP error %s when saving %s to workspace %s' %(e.reason, input['probanno'], input['probanno_workspace'])
                sys.stderr.write('[warning] '+message+'\n')
                self._log(log.WARNING, message)
                time.sleep(15)
        
        # Saving the object failed so raise the last exception that was caught.
        raise e

    def runCalculate(self, input):
        ''' Calculate reaction likelihoods and store them in a RxnProbs typed object.

            A RxnProbs typed object is created in five steps: (1) calculate per-gene
            role likelihoods (2) calculate whole cell role likelihoods (3) calculate
            complex likelihoods (4) calculate reaction likelihoods (5)  save the
            likelihood scores to a RxnProbs typed object.

            @param input: Input parameters for calculate() function
            @return Object info for RxnProbs object
        '''

        # Set log level to DEBUG when verbose parameter is enabled.
        if input['verbose'] and self.logger.get_log_level() < log.DEBUG:
            self.logger.set_log_level(log.DEBUG)

        # Get the static database files.  If the files do not exist and they are downloaded
        # from Shock, it can take a few minutes before they are ready.
        testDataPath = os.path.join(os.environ['KB_TOP'], 'services', 'testdata')
        self.dataParser.getDatabaseFiles(self.logger, testDataPath)

        # Create a workspace client.
        wsClient = Workspace(self.config["workspace_url"], token=self.ctx['token'])

        # Get the ProbAnno object from the specified workspace.
        probannoObjectId = make_object_identity(input["probanno_workspace"], input["probanno"])
        objectList = wsClient.get_objects( [ probannoObjectId ] )
        probannoObject = objectList[0]
        if probannoObject['info'][2] != ProbAnnoType:
            message = "ProbAnno object type %s is not %s for object %s" %(probannoObject['info'][2], ProbAnnoType, probannoObject['info'][1])
            self._log(log.WARNING, message)
            raise WrongVersionError(message)
        genome = probannoObject["data"]["genome"]

        # Create a temporary directory for storing intermediate files when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            workFolder = tempfile.mkdtemp("", "calculate-%s-" %(genome), self.config["work_folder_path"])
            self._log(log.DEBUG2, 'Intermediate files saved in '+workFolder)
        else:
            workFolder = None

        # When a template model is specified, use it to build dictionaries for roles,
        # complexes, and reactions instead of retrieving from static database files.
        complexesToRoles = None
        reactionsToComplexes = None
        if input["template_model"] is not None or input["template_workspace"] is not None:
            if not(input["template_model"] is not None and input["template_workspace"] is not None) :
                message = "Template model workspace is required if template model ID is provided"
                self._log(log.ERR, message)
                raise ValueError(message)

            # Create a dictionary to map a complex to a list of roles and a dictionary
            # to map a reaction to a list of complexes.  The dictionaries are specific to
            # the specified template model instead of covering everything in the central
            # data model.
            complexesToRoles = dict()
            reactionsToComplexes = dict()

            # Get the list of RoleComplexReactions for the template model from the
            # fba modeling service.  The RoleComplexReactions structure has a list
            # of ComplexReactions structures for the given role.  And each ComplexReactions
            # structure has a list of reactions for the given complex.
            fbaClient = fbaModelServices(self.config['fbamodeling_url'], token=self.ctx['token'])
            roleComplexReactionsList = fbaClient.role_to_reactions( { 'templateModel': input['template_model'], 'workspace': input['template_workspace'] } )

            # Build the two dictionaries from the returned list.
            for rcr in roleComplexReactionsList:
                for complex in rcr['complexes']:
                    complexId = re.sub(r'cpx0*(\d+)', r'kb|cpx.\1', complex['name']) # Convert ModelSEED format to KBase format
                    if complexId in complexesToRoles:
                        complexesToRoles[complexId].append(rcr['name'])
                    else:
                        complexesToRoles[complexId] = [ rcr['name'] ]
                    for reaction in complex['reactions']:
                        reactionId = reaction['reaction']
                        if reactionId in reactionsToComplexes:
                            reactionsToComplexes[reactionId].append(complexId)
                        else:
                            reactionsToComplexes[reactionId] = [ complexId ]

        # Calculate per-gene role probabilities.
        roleProbs = self._rolesetProbabilitiesToRoleProbabilities(input, genome, probannoObject["data"]["roleset_probabilities"], workFolder)

        # Calculate whole cell role probabilities.
        # Note - eventually workFolder will be replaced with a rolesToReactions call
        totalRoleProbs = self._totalRoleProbabilities(input, genome, roleProbs, workFolder)

        # Calculate complex probabilities.
        complexProbs = self._complexProbabilities(input, genome, totalRoleProbs, workFolder, complexesToRequiredRoles = complexesToRoles)

        # Calculate reaction probabilities.
        reactionProbs = self._reactionProbabilities(input, genome, complexProbs, workFolder, rxnsToComplexes = reactionsToComplexes)

        # If the reaction probabilities were not calculated using the data from the fba modeling service
        # via the template model, we need to convert from the KBase ID format to the ModelSEED format.
        if input["template_model"] is None:
            reactionList = list()
            for index in range(len(reactionProbs)):
                reactionList.append(reactionProbs[index][0])
            EntityAPI = CDMI_EntityAPI(self.config["cdmi_url"])
            numAttempts = 4
            while numAttempts > 0:
                try:
                    numAttempts -= 1
                    reactionData = EntityAPI.get_entity_Reaction( reactionList, [ "source_id" ] )
                    if len(reactionList) == len(reactionData):
                        numAttempts = 0
                except HTTPError as e:
                    pass
            for index in range(len(reactionProbs)):
                rxnId = reactionProbs[index][0]
                if rxnId in reactionData:
                    reactionProbs[index][0] = reactionData[rxnId]['source_id']
                else:
                    self._log(log.DEBUG2, 'Reaction %s was not found in CDM' %(rxnId))
                    rxnnum = re.sub(r'kb\|rxn.', '', rxnId)
                    reactionProbs[index][0] = 'rxn%05d' %(int(rxnnum))

        # Create a reaction probability object
        objectData = dict()
        objectData["genome"] = probannoObject["data"]["genome"]
        objectData['genome_workspace'] = probannoObject['data']['genome_workspace']
        if input["template_model"] is None:
            objectData['template_model'] = 'None'
        else:
            objectData["template_model"] = input["template_model"]
        if input["template_workspace"] is None:
            objectData['template_workspace'] = 'None'
        else:
            objectData["template_workspace"] = input["template_workspace"]
        objectData["probanno"] = input['probanno']
        objectData['probanno_workspace'] = input['probanno_workspace']
        objectData["id"] = input["rxnprobs"]
        objectData["reaction_probabilities"] = reactionProbs

        objectMetaData = { "num_reaction_probs": len(objectData["reaction_probabilities"]) }
        objectProvData = dict()
        objectProvData['time'] = timestamp(0)
        objectProvData['service'] = os.environ.get('KB_SERVICE_NAME', ServiceName)
        objectProvData['service_ver'] = ServiceVersion
        objectProvData['method'] = 'calculate'
        objectProvData['method_params'] = input.items()
        objectProvData['input_ws_objects'] = [ '%s/%s/%d' %(probannoObject['info'][7], probannoObject['info'][1], probannoObject['info'][4]) ]
        objectSaveData = dict();
        objectSaveData['type'] = RxnProbsType
        objectSaveData['name'] = input["rxnprobs"]
        objectSaveData['data'] = objectData
        objectSaveData['meta'] = objectMetaData
        objectSaveData['provenance'] = [ objectProvData ]
        objectInfo = wsClient.save_objects( { 'workspace': input["rxnprobs_workspace"], 'objects': [ objectSaveData ] } )
        output = objectInfo[0]

        message = 'Saved RxnProbs object %s/%s' %(input['rxnprobs_workspace'], input['rxnprobs'])
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)
        return output

    def _rolesetProbabilitiesToRoleProbabilities(self, input, genome, queryToTuplist, workFolder):
        ''' Compute probability of each role from the rolesets for each query protein.

            At the moment the strategy is to take any set of rolestrings containing
            the same roles and add their probabilities.  So if we have hits to both
            a bifunctional enzyme with R1 and R2, and hits to a monofunctional enzyme
            with only R1, R1 ends up with a greater probability than R2.

            I had tried to normalize to the previous sum but I need to be more careful
            than that (I'll put it on my TODO list) because if you have e.g. one hit
            to R1R2 and one hit to R3 then the probability of R1 and R2 will be unfairly
            brought down due to the normalization scheme.

            @param input: Dictionary of input parameters to calculate() function
            @param genome: Genome ID string
            @param queryToTuplist: Dictionary keyed by query gene of list of tuples with roleset and likelihood
            @param workFolder: Path to directory in which to store temporary files
            @return List of tuples with query gene, role, and likelihood
        '''

        message = 'Started computing role probabilities from roleset probabilities for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        # Start with an empty list.
        roleProbs = list()

        # Iterate over all of the query genes in the dictionary.
        # querygene -> [ (roleset1, likelihood_1), (roleset2, likelihood_2), ...]
        for query in queryToTuplist:
            # This section actually does the conversion of likelihoods.
            # See equation 3 in the paper ("Calculating reaction likelihoods" section).
            queryRolesToProbs = dict()
            for tup in queryToTuplist[query]:
                rolelist = tup[0].split(self.config["separator"])
                # Add up all the instances of each particular role on the list.
                for role in rolelist:
                    if role in queryRolesToProbs:
                        queryRolesToProbs[role] += tup[1]
                    else:
                        queryRolesToProbs[role] = tup[1]

            # Add them to the array.
            for role in queryRolesToProbs:
                roleProbs.append( (query, role, queryRolesToProbs[role]) )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            role_probability_file = os.path.join(workFolder, "%s.roleprobs" %(genome))
            with open(role_probability_file, "w") as handle:
                for tuple in roleProbs:
                    handle.write("%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2]))

        message = 'Finished computing role probabilities from roleset probabilities for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        return roleProbs

    def _totalRoleProbabilities(self, input, genome, roleProbs, workFolder):
        ''' Given the likelihood that each gene has each role, estimate the likelihood
            that the entire ORGANISM has that role.

            To avoid exploding the likelihoods with noise, I just take the maximum
            likelihood of any query gene having a function and use that as the
            likelihood that the function exists in the cell.

            A gene is assigned to a role if it is within DILUTION_PERCENT of the maximum
            probability. DILUTION_PERCENT can be adjusted in the config file. For each
            role the maximum likelihood and the estimated set of genes that perform that
            role are linked with an OR relationship to form a Boolean Gene-Function
            relationship.

            @param ctx Current context object
            @param input Dictionary of input parameters to calculate() function
            @param genome Genome ID string
            @param roleProbs List of tuples with query gene, role, and likelihood
            @param workFolder Path to directory in which to store temporary files
            @return List of tuples with role, likelihood, and estimated set of genes that perform the role
            @raise RoleNotFoundError when role is not placed properly in roleToTotalProb dictionary
        '''

        message = 'Started generating whole-cell role probability file for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        # Find maximum likelihood among all query genes for each role.
        # This is assumed to be the likelihood of that role occurring in the organism as a whole.
        roleToTotalProb = dict()
        for tuple in roleProbs:
            if tuple[1] in roleToTotalProb:
                if float(tuple[2]) > roleToTotalProb[tuple[1]]:
                    roleToTotalProb[tuple[1]] = float(tuple[2])
            else:
                roleToTotalProb[tuple[1]] = float(tuple[2])

        # Get the genes within DILUTION_PERCENT percent of the maximum
        # likelihood and assert that these are the most likely genes responsible for that role.
        # (note - DILUTION_PERCENT is defined in the config file)
        # This produces a dictionary from role to a list of genes
        # See equation 4 in the paper ("Calculating reaction likelihoods" section).
        roleToGeneList = dict()
        for tuple in roleProbs:
            if tuple[1] not in roleToTotalProb:
                message = "Role %s not placed properly in roleToTotalProb dictionary?" %(tuple[1])
                self._log(log.ERR, message)
                raise RoleNotFoundError(message)
            if float(tuple[2]) >= float(self.config["dilution_percent"])/100.0 * roleToTotalProb[tuple[1]]:
                if tuple[1] in roleToGeneList:
                    roleToGeneList[tuple[1]].append(tuple[0])
                else:
                    roleToGeneList[tuple[1]] = [ tuple[0] ]

        # Build the array of total role probabilities.
        totalRoleProbs = list()
        for role in roleToTotalProb:
            gpr = " or ".join(list(set(roleToGeneList[role])))
            # We only need to group these if there is more than one of them (avoids extra parenthesis when computing complexes)
            if len(list(set(roleToGeneList[role]))) > 1:
                gpr = "(" + gpr + ")"
            totalRoleProbs.append( (role, roleToTotalProb[role], gpr ) )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            total_role_probability_file = os.path.join(workFolder, "%s.cellroleprob" %(genome))
            with open(total_role_probability_file, "w") as handle:
                for tuple in totalRoleProbs:
                    handle.write("%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2]))

        message = 'Finished generating whole-cell role probability file for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        return totalRoleProbs

    def _complexProbabilities(self, input, genome, totalRoleProbs, workFolder, complexesToRequiredRoles = None):
        ''' Compute the likelihood of each protein complex from the likelihood of each role.

            A protein complex represents a set functional roles that must all be present
            for a complex to exist.  The likelihood of the existence of a complex is
            computed as the minimum likelihood of the roles within that complex (ignoring
            roles not represented in the subsystems).

            For each protein complex, the likelihood, type, list of roles not in the
            organism, and list of roles not in subsystems is returned.  The type is a
            string with one of the following values:

            CPLX_FULL - All roles found in organism and utilized in the complex
            CPLX_PARTIAL - Only some roles found in organism and only those roles that
                were found were utilized. Note this does not distinguish between not
                there and not represented for roles that were not found
            CPLX_NOTTHERE - Likelihood is 0 because the genes aren't there for any of
                the subunits
            CPLX_NOREPS - Likelihood is 0 because there are no representative genes in
                the subsystems for any of the subunits
            CPLX_NOREPS_AND_NOTTHERE - Likelihood is 0 because some genes aren't there
                for any of the subunits and some genes have no representatives

            @param input: Dictionary of input parameters to calculate() function
            @param genome: Genome ID string
            @param totalRoleProbs: List of tuples with role, likelihood, and estimated set
                of genes that perform the role
            @param workFolder: Path to directory in which to store temporary files
            @param complexesToRequiredRoles: Dictionary keyed by complex ID to the roles
                involved in forming that complex. If it is None we read it from the CDMI
                files we downloaded, otherwise we use the provided dictionary. This is
                included for template model support in the future
            @return List of tuples with complex ID, likelihood, type, list of roles not in
                organism, list of roles not in subsystems, and boolean Gene-Protein
                relationship
        '''

        message = 'Started computing complex probabilities for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        # Get the mapping from complexes to roles if it isn't already provided.
        if complexesToRequiredRoles is None:
            complexesToRequiredRoles, rolesToComplexes = self.dataParser.readComplexRoleFile(self.dataParser.DataFiles['complex_role_file'])
            self._log(log.INFO, 'Found %d complex to role mappings in %s' %(len(complexesToRequiredRoles), self.dataParser.DataFiles['complex_role_file']))

        # Get the subsystem roles (used to distinguish between NOTTHERE and NOREPS).
        otu_fidsToRoles, otu_rolesToFids = self.dataParser.readFidRoleFile(self.dataParser.DataFiles['otu_fid_role_file'])
        allroles = set()
        for fid in otu_fidsToRoles:
            for role in otu_fidsToRoles[fid]:
                allroles.add(role)

        # Build two dictionaries, both keyed by role, one mapping the role to its
        # likelihood and one mapping to the gene list.
        rolesToProbabilities = dict()
        rolesToGeneList = dict()
        for tuple in totalRoleProbs:
            rolesToProbabilities[tuple[0]] = float(tuple[1]) # can skip the float()?
            rolesToGeneList[tuple[0]] = tuple[2]

        # Iterate over complexes and compute complex probabilities from role probabilities.
        # Separate out cases where no genes seem to exist in the organism for the reaction
        # from cases where there is a database deficiency.
        # See equation 5 in the paper ("Calculating reaction likelihoods" section).
        SEPARATOR = self.config["separator"]
        complexProbs = list()
        for cplx in complexesToRequiredRoles:
            allCplxRoles = complexesToRequiredRoles[cplx]
            availRoles = list() # Roles that may have representatives in the query organism
            unavailRoles = list() # Roles that have representatives but that are not apparently in the query organism
            noexistRoles = list() # Roles with no representatives in the subsystems
            for role in complexesToRequiredRoles[cplx]:
                if role not in allroles:
                    noexistRoles.append(role)
                elif role not in rolesToProbabilities:
                    unavailRoles.append(role)
                else:
                    availRoles.append(role)
            TYPE = ""
            GPR = ""
            if len(noexistRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOREPS"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            if len(unavailRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOTTHERE"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            # Some had no representatives and the rest were not found in the cell
            if len(unavailRoles) + len(noexistRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOREPS_AND_NOTTHERE"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            # Otherwise at least one of them is available
            if len(availRoles) == len(allCplxRoles):
                TYPE = "CPLX_FULL"
            elif len(availRoles) < len(allCplxRoles):
                TYPE = "CPLX_PARTIAL_%d_of_%d" %(len(availRoles), len(allCplxRoles))

            # Link individual functions in complex with an AND relationship to form a
            # Boolean Gene-Protein relationship.
#            partialGprList = [ "(" + s + ")" for s in [ rolesToGeneList[f] for f in availRoles ] ]
            partialGprList = [ rolesToGeneList[f] for f in availRoles ]
            GPR = " and ".join( list(set(partialGprList)) )

            if GPR != "" and len(list(set(partialGprList))) > 1:
                GPR = "(" + GPR + ")"

            # Find the minimum probability of the different available roles (ignoring ones
            # that are apparently missing) and call that the complex likelihood.
            minp = 1000
            for role in availRoles:
                if rolesToProbabilities[role] < minp:
                    minp = rolesToProbabilities[role]
            complexProbs.append( (cplx, minp, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            complex_probability_file = os.path.join(workFolder, "%s.complexprob" %(genome))
            with open(complex_probability_file, "w") as handle:
                for tuple in complexProbs:
                    handle.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4], tuple[5]))

        message = 'Finished computing complex probabilities for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)
        return complexProbs

    def _reactionProbabilities(self, input, genome, complexProbs, workFolder, rxnsToComplexes = None):
        ''' Estimate the likelihood of reactions from the likelihood of complexes.

            The reaction likelihood is computed as the maximum likelihood of complexes
            that perform that reaction.

            If the reaction has no complexes it won't even be in this file because of the way
            I set up the call... I could probably change this so that I get a list of ALL reactions
            and make it easier to catch issues with reaction --> complex links in the database.
            Some of the infrastructure is already there (with the TYPE).

            @param input: Dictionary of input parameters to calculate() function
            @param genome: Genome ID string
            @param complexProbs: List of tuples with complex ID, likelihood, type, list of
                roles not in organism, list of roles not in subsystems, and boolean
                Gene-Protein relationship
            @param workFolder: Path to directory in which to store temporary files
            @param rxnsToComplexes: Dictionary keyed by reaction ID to a list of catalyzing
                complexes. If it is None we read it from the CDMI files we downloaded,
                otherwise we use the provided dictionary. This is included for template
                model support in the future
            @return List of tuples with reaction ID, likelihood, reaction type, complex info,
                and gene-protein-reaction relationship
        '''

        message = 'Started computing reaction probabilities for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)

        # Build a dictionary keyed by complex ID of tuples with likelihood, type, and GPR.
        # Note we don't need to use the list of roles not in organism and list of roles
        # not in subsystems.
        # cplx --> {likelihood, type, GPR}
        cplxToTuple = dict()
        for tuple in complexProbs:
            cplxToTuple[tuple[0]] = ( tuple[1], tuple[2], tuple[5] )
        
        # Get the mapping from reactions to complexes if it isn't already provided.
        if rxnsToComplexes is None:
            rxnsToComplexes = self.dataParser.readReactionComplexFile(self.dataParser.DataFiles['reaction_complex_file'])

        # Take the MAXIMUM likelihood of complexes catalyzing a particular reaction
        # and call that the reaction likelihood.
        # See equation 6 in the paper ("Calculating reaction likelihoods" section).
        reactionProbs = list()
        for rxn in rxnsToComplexes:
            TYPE = "NOCOMPLEXES"
            rxnComplexes = rxnsToComplexes[rxn]
            maxProb = 0
            GPR = ""
            complexList = list()
            for cplx in rxnComplexes:
                if cplx in cplxToTuple:
                    # Complex1 (P1; TYPE1) ///Complex2 (P2; TYPE2) ...
                    complexList.append( [ cplx, cplxToTuple[cplx][0], cplxToTuple[cplx][1] ])
                    TYPE = 'HASCOMPLEXES'
            complexString = ''
            if len(complexList) > 0:
                complexList.sort(key=lambda tup: tup[1], reverse=True)
                maxProb = complexList[0][1]
                for complex in complexList:
                    complexString += '%s (%1.4f; %s)%s' %(complex[0], complex[1], complex[2], self.config['separator'])
                complexString = complexString[:-len(self.config['separator'])] # Remove the final separator

            # Iterate separately to get a GPR. We want to apply a cutoff here too to avoid
            # a complex with 80% probability being linked by OR to another with a 5%
            # probability.  For now I've implemented using the same cutoff as we used for
            # which genes go with a role.
            cplxGprs = []
            for cplx in rxnComplexes:
                if cplx in cplxToTuple:
                    if cplxToTuple[cplx][0] < maxProb * float(self.config["dilution_percent"])/100.0:
                        continue
                    cplxGprs.append(cplxToTuple[cplx][2])
            if len(cplxGprs) > 0:
                GPR = " or ".join( list(set(cplxGprs)) )

            # Use a list so that we can modify the reaction IDs if needed to translate to ModelSEED IDs
            reactionProbs.append( [rxn, maxProb, TYPE, complexString, GPR] )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            reaction_probability_file = os.path.join(workFolder, "%s.rxnprobs" %(genome))
            with open(reaction_probability_file, "w") as handle:
                for tuple in reactionProbs:
                    handle.write("%s\t%1.4f\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]))

        message = 'Finished computing reaction probabilities for '+genome
        sys.stdout.write(message+'\n')
        self._log(log.DEBUG, message)
        return reactionProbs

    def _log(self, level, message):
        ''' Log a message to the system log.

            @param level: Message level (INFO, WARNING, etc.)
            @param message: Message text
            @return Nothing
        '''

        # Log the message.
        self.logger.log_message(level, message, self.ctx['client_ip'], self.ctx['user_id'], self.ctx['module'],
                                self.ctx['method'], self.ctx['call_id'])
        return

    def __init__(self, context=None):
        ''' Initialize object.

            @return Nothing
        '''

        # Use the context from the server or build a context when running in NJS.
        # @todo Still not quite sure how context is handled in NJS environment.
        submod = os.environ.get('KB_SERVICE_NAME', ServiceName)
        if context is not None:
            self.ctx = context
        else:
            self.ctx = dict()
            self.ctx['token'] = os.environ['KB_AUTH_TOKEN']
            self.ctx['client_ip'] = '127.0.0.1'
            self.ctx['user_id'] = '-'
            self.ctx['module'] = submod
            self.ctx['method'] = '-'
            self.ctx['call_id'] = '-'

        # Create a logger.
        print os.getenv('KB_DEPLOYMENT_CONFIG')
        self.logger = log.log(submod, ip_address=True, authuser=True, module=True, method=True,
            call_id=True, config=os.getenv('KB_DEPLOYMENT_CONFIG'))
        self.logger.log_message(log.INFO, 'Worker started, version is '+ServiceVersion)

        # Get the configuration.
        self.config = get_config()
        configValues = 'shock_url='+self.config['shock_url']
        configValues += ', userandjobstate_url='+self.config['userandjobstate_url']
        configValues += ', workspace_url='+self.config['workspace_url']
        configValues += ', cdmi_url='+self.config['cdmi_url']
        configValues += ', work_folder_path='+self.config['work_folder_path']
        configValues += ', data_folder_path='+self.config['data_folder_path']
        configValues += ', load_data_option='+self.config['load_data_option']
        configValues += ', separator='+self.config['separator']
        configValues += ', dilution_percent='+self.config['dilution_percent']
        configValues += ', pseudo_count='+self.config['pseudo_count']
        configValues += ', job_queue='+self.config['job_queue']
        configValues += ', search_program='+self.config['search_program']
        configValues += ', search_program_path='+self.config['search_program_path']
        configValues += ', blast_threads='+self.config['blast_threads']
        configValues += ', usearch_accel='+self.config['usearch_accel']
        self.logger.log_message(log.INFO, configValues)

        # Create a DataParser object for working with the static database files.
        self.dataParser = DataParser(self.config)

        return
