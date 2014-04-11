
from biokbase.probabilistic_annotation.Helpers import make_object_identity, make_job_directory, timestamp, ProbAnnoType
from biokbase.probabilistic_annotation.DataParser import checkIfDatabaseFilesExist, DatabaseFiles, readFilteredOtuRoles, parseBlastOutput, readFilteredOtuRoles
from biokbase.userandjobstate.client import UserAndJobState
from biokbase.workspace.client import Workspace
from biokbase import log
from urllib2 import HTTPError
import subprocess
import sys
import os
import shutil
import traceback

# Exception thrown when no features are found in Genome object
class NoFeaturesError(Exception):
    pass

# Exception thrown when blast command failed
class BlastError(Exception):
    pass

# Exception thrown when there are no gene IDs in Genome object
class NoGeneIdsError(Exception):
    pass

''' Worker for long running probabilistic annotation jobs. '''

class ProbabilisticAnnotationWorker:

    def runAnnotate(self, job):

        ''' Run an annotate job to create a ProbAnno typed object.

            A ProbAnno typed object is created in four steps: (1) extract amino acid
            sequences from a Genome typed object to a fasta file, (2) run a BLAST search
            using the amino acid sequences against the subsystem BLAST database,
            (3) calculate annotation likelihood scores for each roleset implied by the
            functions of proteins in subsystems, and (4) save the likelihood scores
            to a ProbAnno typed object.

            The Job dictionary contains three main sections: (1) input parameters to
            the annotate() function, (2) context of server instance running the
            annotate() function, and (3) config variables of server.

            @param job Job dictionary created by server's annotate() function
            @return Nothing (although job is marked as complete)
        '''

        # The input parameters and user context for annotate() were stored in the job data for the job.
        input = job["input"]
        self.ctx = job["context"]
        self.config = job['config']

        status = None

        try:
            # Make sure the database files are available.
            checkIfDatabaseFilesExist(self.config['data_folder_path'])

            # Make sure the job directory exists.
            workFolder = make_job_directory(self.config['work_folder_path'], job['id'])

            # Create a user and job state client and authenticate as the user.
            ujsClient = UserAndJobState(self.config['userandjobstate_url'], token=self.ctx['token'])
    
            # Get the Genome object from the specified workspace.
            try:
                ujsClient.update_job_progress(job['id'], self.ctx['token'], 'getting genome object', 1, timestamp(3600))
            except:
                pass
            wsClient = Workspace(self.config["workspace_url"], token=self.ctx['token'])
            genomeObjectId = make_object_identity(input["genome_workspace"], input["genome"])
            objectList = wsClient.get_objects( [ genomeObjectId ] )
            genomeObject = objectList[0]
            
            # Convert Genome object to fasta file.
            try:
                ujsClient.update_job_progress(job['id'], self.ctx['token'], 'converting Genome object to fasta file', 1, timestamp(3600))
            except:
                pass
            fastaFile = self._genomeToFasta(input, genomeObject, workFolder)
            
            # Run blast using the fasta file.
            try:
                ujsClient.update_job_progress(job['id'], self.ctx['token'], 'running blast', 1, timestamp(3600))
            except:
                pass
            blastResultFile = self._runBlast(input, fastaFile, workFolder)
            
            # Calculate roleset probabilities.
            try:
                ujsClient.update_job_progress(job['id'], self.ctx['token'], 'calculating roleset probabilities', 1, timestamp(300))
            except:
                pass
            rolestringTuples = self._rolesetProbabilitiesMarble(input, blastResultFile, workFolder)
            
            # Build ProbAnno object and store in the specified workspace.
            try:
                ujsClient.update_job_progress(job['id'], self.ctx['token'], 'building ProbAnno object', 1, timestamp(120))
            except:
                pass
            output = self._buildProbAnnoObject(input, genomeObject, blastResultFile, rolestringTuples, workFolder, wsClient)

            # Mark the job as done.
            status = "done"
            tb = None

        except:
            tb = traceback.format_exc()
            sys.stderr.write('\n'+tb)
            status = "failed"
        
        # Mark the job as complete with the given status.
        ujsClient.complete_job(job['id'], self.ctx['token'], status, tb, { })

        # Remove the temporary work directory.
        if not self.config["debug"] or self.config["debug"] == "0":
            try:
                shutil.rmtree(workFolder)
            except OSError:
                # For some reason deleting the directory was failing in production. Rather than have all jobs look like they failed
                # I catch and log the exception here (since the user still gets the same result if the directory remains intact)
                msg = 'Unable to delete temporary directory %s\n' %(workFolder)
                sys.stderr.write('WARNING: '+msg)
                self._log(log.WARN, msg)

        return
        
    def _genomeToFasta(self, input, genomeObject, workFolder):

        ''' Convert a Genome object into an amino-acid FASTA file (for BLAST purposes).

            @param input Dictionary of input parameters to annotate() function
            @param genomeObject Genome typed object from workspace
            @param workFolder Path to directory in which to store temporary files
            @return Path to fasta file with query proteins
        '''

        # Make sure the Genome object has features.
        if "features" not in genomeObject["data"]:
            raise NoFeaturesError("The input Genome object %s/%s has no features. Did you forget to run annotate_genome?\n" %(input["genome_workspace"], input["genome"]))
    
        # Open the fasta file.
        fastaFile = os.path.join(workFolder, "%s.faa" %(input["genome"]))
        fid = open(fastaFile, "w")
        
        # Run the list of features to build the fasta file.
        features = genomeObject["data"]["features"]
        for feature in features:
            # Not a protein-encoding gene
            if "protein_translation" not in feature:
                continue
            myid = feature["id"]
            if "function" in feature:
                function = feature["function"]
            else:
                function = ""
            seq = feature["protein_translation"]
            fid.write(">%s %s\n%s\n" %(myid, function, seq))
        
        fid.close()    
        return fastaFile
        
    def _runBlast(self, input, queryFile, workFolder):

        ''' A simplistic wrapper to BLAST the query proteins against the subsystem proteins.

            @param input Dictionary of input parameters to annotate() function
            @param queryFile Path to fasta file with query proteins
            @param workFolder Path to directory in which to store temporary files
            @return Path to output file from BLAST
        '''

        # Generate path to output file.  Output format 6 is tab-delimited format.
        blastResultFile = os.path.join(workFolder, "%s.blastout" %(input["genome"]))

        args = [ "blastp", "-query", queryFile, 
                 "-db", os.path.join(self.config["data_folder_path"], DatabaseFiles["subsystem_otu_fasta_file"]),
                 "-outfmt", "6", "-evalue", "1E-5",
                 "-num_threads", self.config["blast_threads"],
                 "-out", blastResultFile ]
        cmd = ' '.join(args)
        sys.stderr.write("Started BLAST with command: %s\n" %(cmd))
        status = subprocess.call(args)
        sys.stderr.write("Ended BLAST with command: %s\n" %(cmd))
        if os.WIFEXITED(status):
            if os.WEXITSTATUS(status) != 0:
                raise BlastError("'%s' failed with status %d\n" %(cmd, os.WEXITSTATUS(status)))
        if os.WIFSIGNALED(status):
            raise BlastError("'%s' ended by signal %d\n" %(cmd, os.WTERMSIG(status)))
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
        '''

        sys.stderr.write("Performing marble-picking on rolesets for genome %s..." %(input["genome"]))
    
        # Read in the target roles (this function returns the roles as lists!)
        targetIdToRole, targetRoleToId = readFilteredOtuRoles(self.config)
    
        # Convert the lists of roles into "rolestrings" (sort the list so that order doesn't matter)
        # in order to deal with the case where some of the hits are multi-functional and others only have
        # a single function.
        targetIdToRoleString = dict()
        for target in targetIdToRole:
            stri = self.config["separator"].join(sorted(targetIdToRole[target]))
            targetIdToRoleString[target] = stri

        # Parse the output from BLAST which returns a dictionary keyed by query gene of a list
        # of tuples with target gene and score.
        # query --> [ (target1, score 1), (target 2, score 2), ... ]
        idToTargetList = parseBlastOutput(blastResultFile)
    
        # This is a holder for all of our results which is a dictionary keyed by query gene
        # of a list of tuples with roleset and likelihood.
        # query -> [ (roleset1, likelihood_1), (roleset2, likelihood_2), ...]
        rolestringTuples = dict()

        # For each query gene we calculate the likelihood of each possible rolestring.
        # See equation N in the paper.
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
    #                sys.stderr.write("ERROR: Target ID %s from the BLAST file had no roles in the rolestring dictionary??\n" %(tup[0]))
                    continue
                if rolestring in rolestringToScore:
                    rolestringToScore[rolestring] += (float(tup[1]) ** 2)
                else:
                    rolestringToScore[rolestring] = (float(tup[1]) ** 2)
    
            # Calculate the likelihood that this gene has the given functional annotation.
            # Start with the denominator which is the sum of squares of the log-scores for
            # all possible rolestrings.
            denom = float(self.config["pseudo_count"]) * maxscore
            for stri in rolestringToScore:
                denom += rolestringToScore[stri]

            # The numerators are the sum of squares for each rolestring.
            # Calculate the likelihood for each rolestring and store in the output dictionary.
            for stri in rolestringToScore:
                p = rolestringToScore[stri] / denom
                if query in rolestringTuples:
                    rolestringTuples[query].append( (stri, p) )
                else:
                    rolestringTuples[query] = [ (stri, p) ]
    
        # Save the generated data when debug is turned on.
        if self.config["debug"]:
            rolesetProbabilityFile = os.path.join(workFolder, "%s.rolesetprobs" %(input['genome']))
            fid = open(rolesetProbabilityFile, "w")
            for query in rolestringTuples:
                for tup in rolestringTuples[query]:
                    fid.write("%s\t%s\t%1.4f\n" %(query, tup[0], tup[1]))
            fid.close()
            
        sys.stderr.write("done\n")
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
        '''
    
        sys.stderr.write("Building ProbAnno object %s/%s for genome %s..." %(input["probanno_workspace"], input["probanno"], input["genome"]))

        # Read in the target roles (this function returns the roles as lists!)
        targetToRoles, rolesToTargets = readFilteredOtuRoles(self.config)
        targetToRoleSet = dict()
        for target in targetToRoles:
            stri = self.config["separator"].join(sorted(targetToRoles[target]))
            targetToRoleSet[target] = stri

        # This is a dictionary from query ID to (target, -log E-value) pairs.
        # We just use it to identify whether or not we actually hit anything in the db
        # when searching for the query gene.
        queryToTargetEvals = parseBlastOutput(blastResultFile)
        
        # For each query ID:
        # 1. Identify their rolestring probabilities (these are the first and second elements of the tuple)
        # 2. Iterate over the target genes and identify those with each function (a list of these and their blast scores is
        #    the third element of the tuple) - should be able to set it up to only iterate once over the list.
        # 3. Add that tuple to the JSON file with the key "alternativeFunctions"
    
        # The Genome object data ["features"] is a list of dictionaries. We want to make our data structure and 
        # then add that to the dictionary.  I use the ii in range so I can edit the elements without changes being lost.
    
        objectData = dict()
        objectData["id"] = input["probanno"]
        objectData["genome"] = input["genome"]
        objectData["genome_workspace"] = input["genome_workspace"];
        objectData["roleset_probabilities"] = queryToRolesetProbs;
        objectData["skipped_features"] = []
        
        for ii in range(len(genomeObject["data"]["features"])):
            feature = genomeObject["data"]["features"][ii]
            if "id" not in genomeObject["data"]:
                raise NoGeneIdsError("No gene IDs found in input Genome object %s/%s (this should never happen)" %(input["genome_workspace"], input["genome"]))
            queryid = feature["id"]
    
            # This can happen if I couldn't find hits from that gene to anything in the database. In this case, I'll just skip it.
            # TODO Or should I make an empty object? I should ask Chris.
            if queryid not in queryToRolesetProbs or queryid not in queryToTargetEvals:
                objectData["skipped_features"].append(queryid)
                
        # Store the ProbAnno object in the specified workspace.
        objectMetaData = dict()
        objectMetaData['num_rolesets'] = len(objectData["roleset_probabilities"])
        objectMetaData['num_skipped_features'] = len(objectData["skipped_features"])
        objectProvData = dict()
        objectProvData['time'] = timestamp(0)
        objectProvData['service'] = os.environ['KB_SERVICE_NAME']
        objectProvData['script'] = 'pa-annotate'
        objectProvData['input_ws_objects'] = [ '%s/%s/%d' %(genomeObject['info'][7], genomeObject['info'][1], genomeObject['info'][4]) ]
        objectSaveData = dict()
        objectSaveData['type'] = ProbAnnoType
        objectSaveData['name'] = input["probanno"]
        objectSaveData['data'] = objectData
        objectSaveData['meta'] = objectMetaData
        objectSaveData['provenance'] = [ objectProvData ]
        retryCount = 3
        while retryCount > 0:
            try:
                objectInfo = wsClient.save_objects( { 'workspace': input["probanno_workspace"], 'objects': [ objectSaveData ] } )
                sys.stderr.write("done\n")
                return objectInfo[0]
            except HTTPError as e:
                # Hopefully this is just a temporary glitch, try again since we worked so hard to build the object.
                retryCount -= 1
                self._log(log.WARN, 'HTTP error %s when saving %s to workspace %s' %(e.reason, input['probanno'], input['probanno_workspace']))
        
        # Saving the object failed so raise the last exception that was caught.
        raise e

    ''' Log a message to the system log.

        @param level Message level (INFO, WARNING, etc.)
        @param message Message text
        @return Nothing
    '''

    def _log(self, level, message):
        # Create a logger if this is the first time the method has been called.
        if self.logger is None:
            submod = os.environ.get('KB_SERVICE_NAME', 'ProbabilisticAnnotation')
            self.logger = log.log(submod, ip_address=True, authuser=True, module=True, method=True,
                call_id=True, config=os.getenv('KB_DEPLOYMENT_CONFIG'))

        # Log the message.
        self.logger.log_message(level, message, self.context['client_ip'], self.context['user_id'], self.context['module'],
                                self.context['method'], self.context['call_id'])
        return

    def __init__(self):
        self.logger = None
