module ProbabilisticAnnotation
{

	/* ************************************************************************************* */
	/* SIMPLE ID AND STRING TYPES */
	/* ************************************************************************************* */

	/* Indicates true or false values (false <= 0, true >=1) */
	typedef int bool;
    
	/* A string identifier for a probabilistic annotation object. */
    typedef string probanno_id;
    
	/* A string identifier for a workspace. Any string consisting of alphanumeric characters and "-" is acceptable. */
	typedef string workspace_id;
	
	/* A string indicating the type of an object stored in a workspace. */
	typedef string object_type;
	
	/* A string identifier for a model object. */
	typedef string model_id;
	
	/* ID of an object stored in the workspace. */
	typedef string object_id;
	
	/* Login name of KBase user account to which permissions for workspaces are mapped */
	typedef string username;
	
	/* Exact time for workspace operations. e.g. 2012-12-17T23:24:06 */
	typedef string timestamp;
	
    /* A permanent reference to an object in a workspace. */
    typedef string workspace_ref;
    
    /* Meta data associated with an object stored in a workspace.
	
		object_id id - ID of the object assigned by the user or retreived from the IDserver (e.g. kb|g.0)
		object_type type - type of the object (e.g. Genome)
		timestamp moddate - date when the object was modified by the user (e.g. 2012-12-17T23:24:06)
		int instance - instance of the object, which is equal to the number of times the user has overwritten the object
		timestamp date_created - time at which the alignment was built/loaded in seconds since the epoch
		string command - name of the command last used to modify or create the object
		username lastmodifier - name of the user who last modified the object
		username owner - name of the user who owns (who created) this object
		workspace_id workspace - ID of the workspace in which the object is currently stored
		workspace_ref ref - a 36 character ID that provides permanent undeniable access to this specific instance of this object
		string chsum - checksum of the associated data object
		mapping<string,string> metadata - custom metadata entered for data object during save operation 
	
	*/
	typedef tuple<object_id id, object_type type, timestamp moddate, int instance, string command, username lastmodifier, username owner, workspace_id workspace, workspace_ref ref, string chsum, mapping<string,string> metadata> object_metadata;

    /*********************************************************************************
    Genome object (from CDM) type definition
    *********************************************************************************/

    /* WARNING: This genome object is NOT THE SAME as the one you get when you load a genome
    into the workspaces.
    Rather it is the one you get when you call annotate_genome or cs_to_genome
    */
    typedef string md5;
    typedef list<md5> md5s;
    typedef string genome_id;
    typedef string contig_id;
    typedef string feature_id;
    typedef string feature_type;
    
    /* A region of DNA is maintained as a tuple of four components:

        the contig
        the beginning position (from 1)
        the strand
        the length

        We often speak of "a region".  By "location", we mean a sequence
        of regions from the same genome (perhaps from distinct contigs).
    */
    typedef tuple<contig_id, int begin, string strand,int length> region_of_dna;

    /*
        a "location" refers to a sequence of regions
    */
    typedef list<region_of_dna> location;
    typedef tuple<string comment, string annotator, int annotation_time> annotation;
    typedef tuple<string function, float probability> alt_func;

    typedef structure {
	feature_id id;
	location location;
	feature_type type;
	string function;
	list<alt_func> alternative_functions;
	string protein_translation;
	list<string> aliases;
	list<annotation> annotations;
    } feature;

    typedef structure {
	contig_id id;
	string dna;
    } contig;

    typedef structure {
	genome_id id;
	string scientific_name;
	string domain;
	int genetic_code;
	string source;
	string source_id;
	list<contig> contigs;
	list<feature> features;
    } GenomeObject;
        
	/* ************************************************************************************* */
	/* PROBABILISTIC ANNOTATION DATA TYPES */
	/* ************************************************************************************* */

    /* Annotation probability for an alternative function
    
    	string function - the name of the functional role being annotated to the feature
    	float probability - the probability that the functional role is associated with the feature
	*/    	
    typedef tuple<string function, float probability> FunctionProbability;
    
    /* Alternative functions for each feature
    
        feature_id id - ID of feature the annotation is associated with 
        list<FunctionProbability> alternativeFunctions - list of alternative functions and probabilities
    */
    typedef structure {
		feature_id id;
		list<FunctionProbability> alternativeFunctions;
    } ProbAnnoFeature;
    
    /* Object to carry alternative functions and probabilities for genes in a genome    

        probanno_id id - ID of the probabilistic annotation object    
        genome_id genome - ID of the genome the probabilistic annotation was built for
        workspace_id genome_workspace - ID of the workspace containing genome
        list<ProbAnnoFeature> featureAlternativeFunctions - list of ProbAnnoFeature objects holding alternative functions for features
        mapping<feature_id feature, list<FunctionProbability>> rolesetProbabilities - mapping of features to list of FunctionProbability objects
        list<feature_id> skippedFeatures - list of features in genome with no probability
    */
    typedef structure {
		probanno_id id;
		genome_id genome;
		workspace_id genome_workspace;
		list<ProbAnnoFeature> featureAlternativeFunctions;
		mapping<feature_id feature, list<FunctionProbability> alternativeFunctions> rolesetProbabilities;
		list<feature_id> skippedFeatures;
    } ProbabilisticAnnotation;
    
    /*********************************************** 
                     Function definitions
    ************************************************/

    /* 
	 probanno_workspace: Workspace to which to input the final probabilistic annotation. Required.
	 genome_workspace: Workspace from which the genome object was taken. Required.
	 probanno: ID with which to save the probanno object. Required.
	 GenomeObject: An already-loaded Genome object. Required.
    */
    typedef structure {
	probanno_id probanno;
	GenomeObject genomeObj;
    } annotation_probabilities_input;

    /*
     * Given a genome object populated with genes and annotations, this function adds
     * potential alternative functions with probabilities
     */
    funcdef annotation_probabilities(annotation_probabilities_input) returns (ProbabilisticAnnotation);

    /* Input parameters for the "annotate" function.

       genome_id genome - ID of Genome object
       workspace_id genome_workspace - ID of workspace where Genome object is stored
       probanno_id probanno - ID of ProbAnno object
       workspace_id probanno_workspace - ID workspace where ProbAnno object is saved
       bool overwrite - True to overwrite existing ProbAnno object with same name
       bool debug - True to keep intermediate files for debug purposes
       string auth - Authentication token of KBase user
    */
    typedef structure {
		genome_id genome;
		workspace_id genome_workspace;
		probanno_id probanno;
		workspace_id probanno_workspace;
		bool overwrite;
		bool debug;
		string auth;
    } annotate_params;

    /* This function, rather than using an already-loaded genome object, loads a genome from the specified workspace
       before running the probabilistic annotation algorithm.
     */
    funcdef annotate(annotate_params input) returns (object_metadata output);
    
    /* Input parameters for the "calculate" function.
    
		probanno_id probanno - ID of ProbAnno object
		workspace_id probanno_workspace - ID of workspace where ProbAnno object is stored
		model_id model - ID of Model object
		workspace_id model_workspace - ID of workspace where Model object is saved   
		bool overwrite - True to overwrite existing ProbAnno object with same name
		bool debug - True to keep intermediate files for debug purposes
		string auth - Authentication token of KBase user
    */
    typedef structure {
    	probanno_id probanno;
    	workspace_id probanno_workspace;
    	model_id model;
    	workspace_id model_workspace;
    	bool overwrite;
    	bool debug;
    	string auth;
    } calculate_params;
    
	funcdef calculate(calculate_params input) returns(object_metadata output);
	    	
    /* Input parameters for the "normalize" function.
    
		model_id model - ID of Model object
		workspace_id model_workspace - ID of workspace where Model object is saved   
		bool debug - True to keep intermediate files for debug purposes
		string auth - Authentication token of KBase user
    */
    typedef structure {
    	model_id model;
    	workspace_id model_workspace;
    	bool debug;
    	string auth;
    } normalize_params;
    
	funcdef normalize(normalize_params input) returns(bool success);
	    	
    /* Input parameters for the "generate_data" function.
    
    	string folder - Path to folder for generated data files
    	bool regenerate - True to delete and regenerate existing data files
    	bool delete_only - True to only delete existing data files
    	bool verbose - True to enable verbose output
    */
    typedef structure {
    	string folder;
    	bool regenerate;
    	bool delete_only;
    	bool verbose;
    } generate_data_params;
    
    funcdef generate_data(generate_data_params input) returns(bool success);
};
