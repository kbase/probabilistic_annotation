#BEGIN_HEADER
import sys
import os
import traceback
import json
from biokbase.probabilistic_annotation.DataParser import DataParser, NotReadyError
from biokbase.probabilistic_annotation.Helpers import timestamp, make_object_identity, make_job_directory, ProbAnnoType, RxnProbsType, WrongVersionError, ServiceVersion, ServiceName
from biokbase.probabilistic_annotation.Worker import ProbabilisticAnnotationWorker
from biokbase.workspace.client import Workspace
from biokbase.userandjobstate.client import UserAndJobState
from biokbase import log

#END_HEADER


class ProbabilisticAnnotation:
    '''
    Module Name:
    ProbabilisticAnnotation

    Module Description:
    The purpose of the Probabilistic Annotation service is to provide users with
alternative annotations for genes, each attached to a likelihood score, and to
translate these likelihood scores into likelihood scores for the existence of
reactions in metabolic models.  With the Probabilistic Annotation service:

- Users can quickly assess the quality of an annotation.

- Reaction likelihood computations allow users to estimate the quality of
  metabolic networks generated using the automated reconstruction tools in
  other services.

- Combining reaction likelihoods with gapfilling both directly incorporates
  available genetic evidence into the gapfilling process and provides putative
  gene annotations automatically, reducing the effort needed to search for
  evidence for gapfilled reactions.
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    def _checkInputArguments(self, ctx, input, requiredArgs, defaultArgDict):
        ''' Check that required input arguments are present and set defaults for non-required arguments.
        
            If a key in defaultArgDict is not found in the input it is added with the
            specified default value.

            @param ctx Current context object
            @param input Dictionary keyed by argument name of argument values
            @param requiredArgs List of required arguments in the input dictionary
            @param defaultArgDict Dictionary keyed by argument name of default argument values
            @return Dictionary keyed by argument name of argument values (including default arguments)
            @raise ValueError when required argument is not found
        '''
        if requiredArgs is not None:
            for arg in requiredArgs:
                if arg not in input:
                    message = "Required argument %s not found" %(arg)
                    ctx.log_err(message)
                    raise ValueError(message)
        if defaultArgDict is not None:
            for arg in defaultArgDict:
                if arg in input:
                    continue
                else:
                    input[arg] = defaultArgDict[arg]

        return input

    def _checkDatabaseFiles(self, ctx):
        ''' Check the status of the static database files.

            @param ctx Current context object
            @return Nothing
            @raise NotReadyError if the database has not been loaded correctly.
        '''
        try:
            status = self.dataParser.readStatusFile()
            if status != "ready":
                message = "Static database files are not ready.  Current status is '%s'." %(status)
                ctx.log_err(message)
                raise NotReadyError(message)
        except IOError:
            message = "Static database files are not ready.  Failed to open status file '%s'." %(self.dataParser.StatusFiles['status_file'])
            ctx.log_err(message)
            raise NotReadyError(message)
        return

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        ''' Constructor for ProbabilisticAnnotation object.

            @param config Dictionary of configuration variables
            @return Nothing
            @raise ValueError when a valid configuration was not provided
        '''
        if config == None:
            # There needs to be a config for the server to work.
            raise ValueError("__init__: A valid configuration was not provided.  Check KB_DEPLOYMENT_CONFIG and KB_SERVICE_NAME environment variables.")
        else:
            self.config = config
        
        submod = os.environ.get('KB_SERVICE_NAME', ServiceName)
        self.mylog = log.log(submod, ip_address=True, authuser=True, module=True, method=True,
            call_id=True, config=os.getenv('KB_DEPLOYMENT_CONFIG'))
        self.mylog.log_message(log.NOTICE, 'Server started, version is '+ServiceVersion)
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
        self.mylog.log_message(log.NOTICE, configValues)

        # Create a DataParser object for working with the static database files (the
        # data folder is created if it does not exist).
        self.dataParser = DataParser(self.config)

        # Get the static database files.  If the files do not exist and they are downloaded
        # from Shock, it can take a few minutes before the server is ready.
        testDataPath = os.path.join(os.environ['KB_SERVICE_DIR'], 'testdata')
        self.config['load_data_option'] = self.dataParser.getDatabaseFiles(self.mylog, testDataPath)

        # Validate the value of the job_queue variable.  Currently the only supported value is 'local'.
        # Force it to a valid value to avoid an error trying to submit a job later.
        if self.config['job_queue'] != 'local':
            self.mylog.log_message(log.NOTICE, 'Configuration variable job_queue='+self.config['job_queue']+' switched to local')
            self.config['job_queue'] = 'local'
        #END_CONSTRUCTOR
        pass

    def version(self, ctx):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN version
        ''' Return the name and version number of the service.

            @param ctx Current context object
            @return List with service name string and version number string
        '''

        returnVal = [ os.environ.get('KB_SERVICE_NAME'), ServiceVersion ]
        #END version

        # At some point might do deeper type checking...
        if not isinstance(returnVal, list):
            raise ValueError('Method version return value ' +
                             'returnVal is not type list as required.')
        # return the results
        return [returnVal]

    def annotate(self, ctx, input):
        # ctx is the context object
        # return variables are: jobid
        #BEGIN annotate
        ''' Compute probabilistic annotations from the specified genome object.

            The input dictionary must contain the following keys:
            genome: Name of genome object
            genome_workspace: Workspace from which to grab the Genome object
            probanno: Name of probanno object to output
            probanno_workspace: Workspace to which to save the ProbAnno object

            The following keys are optional:
            verbose: Print lots of messages on the progress of the algorithm

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return Job ID of job started to compute annotation likelihoods
            @raise NotReadyError when static database files are not ready
        '''

        input = self._checkInputArguments(ctx, input, 
                                          [ "genome", "genome_workspace", "probanno", "probanno_workspace"],
                                          { "verbose" : False }
                                          )
        
        # Make sure the static database files are ready.
        self._checkDatabaseFiles(ctx)
        
        # Set log level to DEBUG when verbose parameter is enabled.
        if input['verbose']:
            ctx.set_log_level(log.DEBUG)

        # Make sure the Genome object is available.
        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        genomeIdentity = make_object_identity(input['genome_workspace'], input['genome'])
        wsClient.get_object_info( [ genomeIdentity ], 0 )

        # Create a user and job state client and authenticate as the user.
        ujsClient = UserAndJobState(self.config['userandjobstate_url'], token=ctx['token'])

        # Create a job to track running probabilistic annotation.
        description = 'pa-annotate for genome %s to probanno %s for user %s' %(input['genome'], input['probanno'], ctx['user_id'])
        progress = { 'ptype': 'task', 'max': 5 }
        jobid = ujsClient.create_and_start_job(ctx['token'], 'initializing', description, progress, timestamp(3600))
        ctx.log_info('Job '+jobid+' started for genome '+input['genome']+' to probanno '+input['probanno'])

        try:
            # Run the job on the local machine.
            # Note that the constructor validated the config variable and 'local' is the only valid value.
            if self.config['job_queue'] == 'local':
                # Create working directory for job and build file names.
                jobDirectory = make_job_directory(self.config['work_folder_path'], jobid)
                jobDataFilename = os.path.join(jobDirectory, 'jobdata.json')
                outputFilename = os.path.join(jobDirectory, 'stdout.log')
                errorFilename = os.path.join(jobDirectory, 'stderr.log')

                # Save data required for running the job.
                jobData = { 'id': jobid, 'input': input, 'context': ctx, 'config': self.config }
                json.dump(jobData, open(jobDataFilename, "w"), indent=4)

                # Start worker to run the job.
                jobScript = os.path.join(os.environ['KB_TOP'], 'bin/pa-runjob')
                cmdline = "nohup %s %s >%s 2>%s &" %(jobScript, jobDirectory, outputFilename, errorFilename)
                status = os.system(cmdline)
                ctx.log_info('Job %s is running on local host, status %d' %(jobid, status))

        except Exception as e:
            # Mark the job as failed.
            tb = traceback.format_exc()
            ctx.log_err(tb)
            ujsClient.complete_job(jobid, ctx['token'], 'failed', tb, { })
        #END annotate

        # At some point might do deeper type checking...
        if not isinstance(jobid, basestring):
            raise ValueError('Method annotate return value ' +
                             'jobid is not type basestring as required.')
        # return the results
        return [jobid]

    def calculate(self, ctx, input):
        # ctx is the context object
        # return variables are: output
        #BEGIN calculate
        ''' Calculate reaction probabilities from a probabilistic annotation.

            The input dictionary must contain the following keys:
            probanno: Name of ProbAnno object to input
            probanno_workspace: Workspace from which to grab the ProbAnno object
            rxnprobs: Name of RxnProbs object
            rxnprobs_workspace: Workspace to which to save the RxnProbs object

            The following keys are optional:
            verbose: Print lots of messages on the progress of the algorithm
            template_model: Name of TemplateModel object
            template_workspace: Workspace from which to grab TemplateModel object

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return Object info for RxnProbs object
            @raise ValueError when template_workspace input argument is not specified
            @raise NotReadyError when static database files are not ready
        '''

        # Sanity check on input arguments
        input = self._checkInputArguments(ctx, input, 
                                          ["probanno", "probanno_workspace", "rxnprobs", "rxnprobs_workspace"], 
                                          { "verbose" : False ,
                                            "template_model" : None,
                                            "template_workspace" : None
                                          }
                                         )

        # Make sure the static database files are ready.
        self._checkDatabaseFiles(ctx)
        
        # Create a Worker object and calculate reaction probabilities.
        worker = ProbabilisticAnnotationWorker(ctx)
        output = worker.runCalculate(input)
        #END calculate

        # At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method calculate return value ' +
                             'output is not type list as required.')
        # return the results
        return [output]

    def get_rxnprobs(self, ctx, input):
        # ctx is the context object
        # return variables are: output
        #BEGIN get_rxnprobs
        ''' Convert a reaction probability object into a human-readable table.

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return List of reaction_probability tuples
            @raise WrongVersionError when RxnProbs object version number is invalid
        '''

        # Sanity check on input arguments
        input = self._checkInputArguments(ctx, input, 
                                          [ "rxnprobs", "rxnprobs_workspace" ], 
                                          { 'rxnprobs_version': None, 'sort_field': 'rxnid' }
                                          )

        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        rxnProbsObjectId = make_object_identity(input["rxnprobs_workspace"], input["rxnprobs"], input['rxnprobs_version'])
        objectList = wsClient.get_objects( [ rxnProbsObjectId ] )
        rxnProbsObject = objectList[0]
        if rxnProbsObject['info'][2] != RxnProbsType:
            message = 'RxnProbs object type %s is not %s for object %s' %(rxnProbsObject['info'][2], RxnProbsType, rxnProbsObject['info'][1])
            ctx.log_err(message)
            raise WrongVersionError(message)
        output = rxnProbsObject["data"]["reaction_probabilities"]
        if input['sort_field'] == 'rxnid':
            output.sort(key=lambda tup: tup[0])
        elif input['sort_field'] == 'probability':
            output.sort(key=lambda tup: tup[1], reverse=True)
        #END get_rxnprobs

        # At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method get_rxnprobs return value ' +
                             'output is not type list as required.')
        # return the results
        return [output]

    def get_probanno(self, ctx, input):
        # ctx is the context object
        # return variables are: output
        #BEGIN get_probanno
        ''' Convert a probabilistic annotation object into a human-readbable table.

            @param ctx Current context object
            @param input Dictionary with input parameters for function
            @return Dictionary keyed by gene to a list of tuples with roleset and likelihood
            @raise WrongVersionError when ProbAnno object version number is invalid
        '''

        input = self._checkInputArguments(ctx, input,
                                          ['probanno', 'probanno_workspace'],
                                          { 'probanno_version': None }
                                          )

        wsClient = Workspace(self.config["workspace_url"], token=ctx['token'])
        probAnnoObjectId = make_object_identity(input["probanno_workspace"], input["probanno"], input['probanno_version'])
        objectList = wsClient.get_objects( [ probAnnoObjectId ] )
        probAnnoObject = objectList[0]
        if probAnnoObject['info'][2] != ProbAnnoType:
            message = 'ProbAnno object type %s is not %s for object %s' %(probAnnoObject['info'][2], ProbAnnoType, probAnnoObject['info'][1])
            ctx.log_err(message)
            raise WrongVersionError(message)
        output = probAnnoObject["data"]["roleset_probabilities"]

        #END get_probanno

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method get_probanno return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
