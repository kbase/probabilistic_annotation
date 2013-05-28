module ProbabilisticAnnotation
{

	/* ************************************************************************************* */
	/* SIMPLE ID AND STRING TYPES */
	/* ************************************************************************************* */

	/* Indicates true or false values (false <= 0, true >=1) */
	typedef int bool;
    
	/* A string identifier for a probabilistic annotation object. */
    typedef string probanno_id;
    
    /* A string identifier for a template model object. */
    typedef string template_id;

	/* A string identifier for a reaction probabilities object. */
	typedef string rxnprobs_id;
	
	/* A string identifier for a genome. */    
    typedef string genome_id;
    
    /* A string identifier for a feature. */
    typedef string feature_id;
    
	/* A string identifier for a workspace. Any string consisting of alphanumeric characters and "-" is acceptable. */
	typedef string workspace_id;
	
	/* A string indicating the type of an object stored in a workspace. */
	typedef string object_type;
	
	/* A string identifier for a reaction object. */
	typedef string reaction_id;
	
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

	/* ************************************************************************************* */
	/* PROBABILISTIC ANNOTATION DATA TYPES */
	/* ************************************************************************************* */

    /* Annotation probability for an alternative function
    
    	string function - the name of the functional role being annotated to the feature
    	float probability - the probability that the functional role is associated with the feature
    	string functionMD5 - hash to let us know if anything has changed
	*/    	
    typedef tuple<string function, float probability, string functionMD5> FunctionProbability;
    
    /* Object to carry alternative functions and probabilities for genes in a genome    

        probanno_id id - ID of the probabilistic annotation object    
        genome_id genome - ID of the genome the probabilistic annotation was built for
        workspace_id genome_workspace - ID of the workspace containing genome
        mapping<feature_id, list<FunctionProbability>> rolesetProbabilities - mapping of features to list of alternative FunctionProbability objects
        list<feature_id> skippedFeatures - list of features in genome with no probability
    */
    typedef structure {
		probanno_id id;
		genome_id genome;
		workspace_id genome_workspace;
		mapping<feature_id, list<FunctionProbability>> rolesetProbabilities;
		list<feature_id> skippedFeatures;
    } ProbabilisticAnnotation;
    
    /* Data structure to hold probability of a reaction
    
    	reaction_id reaction - ID of the reaction
    	float probability - Probability of the reaction
    	string type - Type of complexes ("HASCOMPLEXES" or "NOCOMPLEXES")
    	string complex_info - Detailed information on complexes
    	string gene_list - List of genes most likely to be attached to reaction
    	
    */
    typedef tuple<reaction_id reaction, float probability, string type, string complex_info, string gene_list> ReactionProbability;
    
    /* Object to hold reaction probabilities for a genome.
    
    	genome_id genome;
    	list<ReactionProbability> reactionProbabilities;
    	
    */
    typedef structure {
    	genome_id genome;
    	list<ReactionProbability> reactionProbabilities;
    } RxnProbs;

	/* ************************************************************************************* */
	/* FUNCTION DEFINITIONS */
	/* ************************************************************************************* */

    /* Input parameters for the "annotate" function.

       genome_id genome - ID of Genome object
       workspace_id genome_workspace - ID of workspace where Genome object is stored
       probanno_id probanno - ID of ProbAnno object
       workspace_id probanno_workspace - ID workspace where ProbAnno object is saved
       bool overwrite - True to overwrite existing ProbAnno object with same name
       bool debug - True to keep intermediate files for debug purposes
	   bool verbose - True to print verbose messages
       string auth - Authentication token of KBase user
    */
    typedef structure {
		genome_id genome;
		workspace_id genome_workspace;
		probanno_id probanno;
		workspace_id probanno_workspace;
		bool overwrite;
		bool debug;
		bool verbose;
		string auth;
    } annotate_params;

    funcdef annotate(annotate_params input) returns (string jobid);
    
    /* Input parameters for the "calculate" function.
    
		probanno_id probanno - ID of ProbAnno object
		workspace_id probanno_workspace - ID of workspace where ProbAnno object is stored
		bool debug - True to keep intermediate files for debug purposes
		bool verbose - True to print verbose messages
		string auth - Authentication token of KBase user
    */
    typedef structure {
    	probanno_id probanno;
    	workspace_id probanno_workspace;
		template_id template_model;
		workspace_id template_model_workspace;
		rxnprobs_id rxnprobs;
		workspace_id rxnprobs_workspace;
    	bool debug;
    	bool verbose;
    	string auth;
    } calculate_params;
    
    /* Compute reaction probabilities from probabilistic annotation and a template model.
    Returns metadata for the reaction probability object 
    */
    funcdef calculate(calculate_params input) returns(object_metadata output);
	    	
};
