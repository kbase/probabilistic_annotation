module ProbabilisticAnnotation
{
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
        
    /* Data structures to hold a single annotation probability for a single gene
    
    feature_id feature - feature the annotation is associated with
    string function - the name of the functional role being annotated to the feature
    float probability - the probability that the functional role is associated with the feature

    */
    typedef tuple<feature_id feature, string function,float probability> annotationProbability; 
    
    typedef string probanno_id;
    typedef tuple<string function, float probability> alt_func;
    
    /*
        Object to carry alternative functions for each feature
    
        feature_id id
        ID of the feature. Required.
    
        string function
        Primary annotated function of the feature in the genome annotation. Required.
    
        list<alt_func> alternative_functions
        List of tuples containing alternative functions and probabilities. Required.
    */
    typedef structure {
	feature_id id;
	list<alt_func> alternative_functions;
    } ProbAnnoFeature;
    
    /* Object to carry alternative functions and probabilities for genes in a genome    

        probanno_id id - ID of the probabilistic annotation object. Required.    
        genome_id genome - ID of the genome the probabilistic annotation was built for. Required.
        list<ProbAnnoFeature> featureAlternativeFunctions - List of ProbAnnoFeature objects holding alternative functions for features. Required.    
	workspace - ID of the workspace from which the genome ID was taken. Required.
    */

    typedef string workspace_id;

    typedef structure {
	probanno_id id;
	genome_id genome;
	list<ProbAnnoFeature> featureAlternativeFunctions;
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

    /* I will implement getting the genome from the workspace later. For now its a headache that is too much
       beyond all the other stuff I have to change here.

       This genome_id is a CDM genome_id.
    */
    typedef structure {
	workspace_id genome_workspace;
	probanno_id probanno;
	genome_id genome;
    } annotation_probabilities_ids_input;

    /* This function, rather than using an already-loaded genome object, loads a genome from the specified workspace
       before running the probabilistic annotation algorithm.
     */
    funcdef annotation_probabilities_id(annotation_probabilities_ids_input) returns (ProbabilisticAnnotation);
};
