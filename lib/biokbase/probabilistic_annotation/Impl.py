#BEGIN_HEADER
from biokbase.probabilistic_annotation.MyImpl import *
#END_HEADER

'''

Module Name:
ProbabilisticAnnotation

Module Description:


'''
class ProbabilisticAnnotation:

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    def __init__(self, config): #config contains contents of config file in hash or 
                                #None if it couldn't be found
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass

    def annotation_probabilities(self, annotation_probabilities_input):
        # self.ctx is set by the wsgi application class
        # return variables are: returnVal
        #BEGIN annotation_probabilities
        #END annotation_probabilities

        #At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method annotation_probabilities return value returnVal is not type dict as required.')
        # return the results
        return [ returnVal ]
        
    def annotate(self, input):
        # self.ctx is set by the wsgi application class
        # return variables are: output
        #BEGIN annotation_probabilities_id
        
        # TODO Get this value from configuration variable for data directory
        dataFolder = "data"
        
        # Create a workspace client.
        wsClient = workspaceService(WORKSPACE_URL)
        
        # Get the Genome object from the specified workspace.
        getObjectParams = { "type": "Genome", "id": input.genome, "workspace": input.genome_workspace, "auth": input.auth }
        genomeObject = wsClient.get_object(getObjectParams)
        
        # Create a temporary directory for storing blast input and output files.
        workFolder = tempfile.mkdtemp("", "%s-" %(input.genome), dataFolder)
        
        # Convert Genome object to fasta file.
        fastaFile = genomeToFasta(input, genomeObject, workFolder)
        
        # Run blast using the fasta file.
        blastResultFile = runBlast(input, fastaFile, workFolder, dataFolder)
        
        # Calculate roleset probabilities.
        rolestringTuples = rolesetProbabilitiesMarble(input, blastResultFile, workFolder, dataFolder)
        
        # Build ProbAnno object and store in the specified workspace.
        output = buildProbAnnoObject(input, genomeObject, blastResultFile, rolestringTuples, workFolder, dataFolder, wsClient)
        
        # Remove the temporary directory.
        if input.debug == False:
            shutil.rmtree(workFolder)
            
        return output
        #END annotation_probabilities_id

        #At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method annotation_probabilities_id return value output is not type list as required.')
        # return the results
        return [ output ]
        
    def calculate(self, input):
        # self.ctx is set by the wsgi application class
        # return variables are: output
        #BEGIN calculate
        
        # TODO Get this value from configuration variable for data directory
        dataFolder = "data"
        
        # Create a workspace client.
        wsClient = workspaceService(WORKSPACE_URL)
        
        # Get the ProbAnno object from the specified workspace.
        getObjectParams = { "type": "ProbAnno", "id": input.probanno, "workspace": input.probanno_workspace, "auth": input.auth }
        probannoObject = wsClient.get_object(getObjectParams)
        input.genome = probannoObject["data"]["genome"]
        input.genome_workspace = probannoObject["data"]["genome_workspace"]    
        sys.stderr.write("%s/%s\n" %(input.genome_workspace, input.genome))
        
        # Create a temporary directory for storing intermediate files. Only used when debug flag is on.
        workFolder = tempfile.mkdtemp("", "%s-" %(input.genome), dataFolder)
        
        # Calculate per-gene role probabilities.
        roleProbs = rolesetProbabilitiesToRoleProbabilities(input, probannoObject["data"]["rolesetProbabilities"], workFolder)
    
        # Calculate whole cell role probabilities.
        totalRoleProbs = totalRoleProbabilities(input, roleProbs, workFolder)
        
        # Calculate complex probabilities.
        complexProbs = complexProbabilities(input, totalRoleProbs, workFolder, dataFolder)
        
        # Calculate reaction probabilities.
        reactionProbs = reactionProbabilities(input, complexProbs, workFolder, dataFolder)
    
        # Create the model with reaction probabilities.
        output = buildModelObject(input, reactionProbs)
    
        return output
        #END calculate

        #At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError('Method calculate return value output is not type list as required.')
        # return the results
        return [ output ]
        
    def normalize(self, input):
        # self.ctx is set by the wsgi application class
        # return variables are: success
        #BEGIN normalize
        
        ### MBM need to clean up input parameters and wrapper in MyImpl.
        sys.stderr.write("ready to normalize!\n")
        # Create a workspace client.
        wsClient = workspaceService(WORKSPACE_URL)
        
        # Get the Model object from the specified workspace.
        getObjectParams = { "type": "Model", "id": input.model, "workspace": input.model_workspace, "auth": input.auth }
        modelObject = wsClient.get_object(getObjectParams)
        
        # Calculate the metabolite weights and update the Model object.
        success = metaboliteWeights(input, modelObject["data"])
        
        # Save the updated Model object.
        if success:
            saveObjectParams = { "type": "Model", "id": input.model, "workspace": input.model_workspace,
                                 "data": modelObject["data"], 
                                 "command": "pa-normalize",  "auth": input.auth }
            metadata = wsClient.save_object(saveObjectParams)
        
        #END normalize

        #At some point might do deeper type checking...
        if not isinstance(success, int):
            raise ValueError('Method normalize return value success is not type int as required.')
        # return the results
        return [ success ]
        
    def generate_data(self, input):
        # self.ctx is set by the wsgi application class
        # return variables are: success
        #BEGIN generate_data
        # When regenerating or deleting the data, remove all of the files first.
        if input.regenerate or input.delete_only:
            safeRemove(OTU_ID_FILE, input.folder)
            safeRemove(SUBSYSTEM_FID_FILE, input.folder)
            safeRemove(DLIT_FID_FILE, input.folder)
            safeRemove(CONCATINATED_FID_FILE, input.folder)
            safeRemove(SUBSYSTEM_OTU_FID_ROLES_FILE, input.folder)
            safeRemove(SUBSYSTEM_OTU_FASTA_FILE, input.folder)
            safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".psq", input.folder) 
            safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".pin", input.folder)
            safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".phr", input.folder)
            safeRemove(COMPLEXES_ROLES_FILE, input.folder)
            safeRemove(REACTION_COMPLEXES_FILE, input.folder)
        
        # Our job is done if all we want to do is delete files.
        if input.delete_only:
            return(1)
        
        sys.stderr.write("Generating requested data:....\n")
        
        # Get lists of OTUs
        sys.stderr.write("OTU data...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            otus, prokotus = readOtuData(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            otus, prokotus = getOtuGenomeIds(MINN, COUNT)
            writeOtuData(otus, prokotus, input.folder)
        sys.stderr.write("done\n")
        
        # Get a list of subsystem FIDs
        sys.stderr.write("List of subsystem FIDS...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            sub_fids = readSubsystemFids(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            sub_fids = subsystemFids(MINN, COUNT)
            writeSubsystemFids(sub_fids, input.folder)
        sys.stderr.write("done\n")
        
        # Get a list of Dlit FIDSs
        # We include these because having them greatly expands the
        # number of roles for which we have representatives.
        sys.stderr.write("Getting a list of DLit FIDs...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            dlit_fids = readDlitFids(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            dlit_fids = getDlitFids(MINN, COUNT)
            writeDlitFids(dlit_fids, input.folder)
        sys.stderr.write("done\n")
        
        # Concatenate the two FID lists before filtering
        # (Note - doing so after would be possible as well but
        # can lead to the same kinds of biases as not filtering
        # the subsystems... I'm not sure the problem would
        # be as bad for these though)
        sys.stderr.write("Combining lists of subsystem and DLit FIDS...")
        fn = os.path.join(input.folder, CONCATINATED_FID_FILE)
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            all_fids = set()
            for line in open(fn, "r"):
                all_fids.add(line.strip("\r\n"))
            all_fids = list(all_fids)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            all_fids = list(set(sub_fids + dlit_fids))
            f = open(fn, "w")
            for fid in all_fids:
                f.write("%s\n" %(fid))
            f.close()
        sys.stderr.write("done\n")
        
        # Identify roles for the OTU genes
        sys.stderr.write("Roles for un-filtered list...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            all_fidsToRoles, all_rolesToFids = readAllFidRoles(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            all_fidsToRoles, all_rolesToFids = fidsToRoles(all_fids)
            writeAllFidRoles(all_fidsToRoles, input.folder)
        sys.stderr.write("done\n")
        
        # Filter the subsystem FIDs by organism. We only want OTU genes.
        # Unlike the neighborhood analysis, we don't want to include only 
        # prokaryotes here.
        sys.stderr.write("Filtered list by OTUs...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            otu_fidsToRoles, otu_rolesToFids = readFilteredOtuRoles(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            otudict = getOtuGenomeDictionary(MINN, COUNT)
            otu_fidsToRoles, otu_rolesToFids, missing_roles = filterFidsByOtusBetter(all_fidsToRoles, all_rolesToFids, otudict)
            writeFilteredOtuRoles(otu_fidsToRoles, input.folder)
        sys.stderr.write("done\n")
        
        # Generate a FASTA file for the fids in fidsToRoles
        sys.stderr.write("Subsystem FASTA file...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            readSubsystemFasta(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys())
            writeSubsystemFasta(fidsToSeqs, input.folder)
        sys.stderr.write("done\n")
        
        # Complexes --> Roles
        # Needed to go from annotation likelihoods
        # to reaction likelihoods
        # Note that it is easier to go in this direction 
        #    Because we need all the roles in a complex to get the probability of that complex.
        sys.stderr.write("Complexes to roles...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            complexToRequiredRoles = readComplexRoles(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(MINN, COUNT)
            writeComplexRoles(complexToRequiredRoles, input.folder)
        sys.stderr.write("done\n")
        
        # reaction --> complex
        # Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.    
        sys.stderr.write("Reactions to complexes...")
        try:
            if input.verbose:
                sys.stderr.write("reading from file...")
            rxnToComplexes = readReactionComplex(input.folder)
        except IOError:
            if input.verbose:
                sys.stderr.write("failed...generating file...")
            rxnToComplexes, complexesToReactions = reactionComplexLinks(MINN, COUNT)
            writeReactionComplex(rxnToComplexes, input.folder)
        sys.stderr.write("done\n")
        
        sys.stderr.write("Data gathering done...\n")
        return(1)
        #END generate_data

        #At some point might do deeper type checking...
        if not isinstance(success, int):
            raise ValueError('Method generate_data return value success is not type int as required.')
        # return the results
        return [ success ]
        
