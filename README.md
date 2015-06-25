Probabilistic Annotation
========================

The purpose of the Probabilistic Annotation service is to provide users with
alternative annotations for genes, each attached to a likelihood score.  The
annotation likelihood scores are then translated into likelihood scores for the
existence of reactions in metabolic models.  With the Probabilistic Annotation
service:

* Users can quickly assess the quality of an annotation.

* Reaction likelihoods estimate the quality of metabolic networks generated
  using the automated reconstruction tools in other services.

* Combining reaction likelihoods with gap filling both directly incorporates
  available genetic evidence into the gap filling process and provides putative
  gene annotations automatically, reducing the effort needed to search for
  evidence for gap filled reactions.

Details on the algorithm are available in [Likelihood-Based Gene Annotations
for Gap Filling and Quality Assessment in Genome-Scale Metabolic
Models](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003882)

Configuration Variables
-----------------------

The Probabilistic Annotation server supports the following configuration variables:

* **default_url**: The default URL for this server's service endpoint.
* **cdmi_url**: URL for the central data model interface service endpoint.  The CDMI
  service is used to build the static database files.
* **shock_url**: URL for the Shock service endpoint.  The Shock service is used to store
  the static database files in a common location.
* **workspace_url**: URL for the Workspace service endpoint.  The Workspace service is
  used to store typed objects created by the service.
* **fbamodeling_url**: URL for the FBA Modeling service endpoint.  The FBA Modeling
  service is used to map roles to reactions when using model templates.
* **userandjobstate_url**: URL of the User and Job State service endpoint.  The UJS
  service is used to manage the jobs that create ProbAnno typed objects.
* **kegg_url**: URL of the Kyoto Encyclopedia of Genes and Genomes API service endpoint.
  The KEGG API service is used update local reaction and enzyme databases and to
  download amino acid sequences for genes when building the static database files.
* **work_folder_path**: Path to work folder containing sub-folders for running jobs.
  The intermediate files created by the annotate() and calculate() methods
  are stored in this location.
* **data_folder_path**: Path to data folder containing static database files.  A server
  can use different static database files by using different folders.
* **load_data_option**: Control how static database files are handled when starting
  service. Valid values are "shock" to load static files from Shock or "preload"
  to use static database files preloaded in the data folder.
* **sources**: List of sources for generating static database files.  The supported
  sources are "cdm" for the KBase central data model and "kegg" for the Kyoto
  Encyclopedia of Genes and Genomes.  Sources are separated by a colon in the list.
  Default value is "cdm".
* **separator**: Character string not found in any roles and used as separator between
  elements in lists. Default value is '///'.
* **dilution_percent**: Percentage of the maximum likelihood to use as a threshold
  to consider other genes as having a particular function aside from the one with
  greatest likelihood. Default value is 80.
* **pseudo_count**: Value used to dilute the likelihoods of annotations for annotations
  with weak homology to the query. Default value is 40.
* **search_program**: Search program for getting log scores of query genes in organism
  against all genes in high-confidence gene annotation database. Valid values are
  "blastp" or "usearch".  Default value is "blastp".
* **search_program_path**: Path to search program.  Use a fully-qualified path name.
* **blast_threads**: Number of threads to use when running search program.
* **search_program_evalue**: Value to use for the search program -evalue parameter.
  Default value is "1E-5".
* **usearch_accel**: Value to use for the -accel parameter of usearch program.  The value
  is a number between 0 and 1 that tunes search speed against sensitivity. Default
  value is 0.33.

Static Database Files
---------------------

The static database files contain pre-processed data that is used by the Probabilistic
Annotation server for building ProbAnno and RxnProbs typed objects.

The static database files come from multiple sources.  The default source is the
KBase central data model and includes these intermediate files:

* **CDM_OTU_GENOME_IDS**: List of representative OTU genome IDs. Each line has one
  field that is a genome ID in KBase format that identifies a genome marked as the
  representative of its OTU.
* **CDM_SUBSYSTEM_FID**: List of feature IDs from SEED subsystems. Each line has one
  field that is the feature ID in KBase format.
* **CDM_DLIT_FID**: List of feature IDs with direct literature support.  Each line
  has one field that is the feature ID in KBase format.
* **CDM_ALL_FID**: Concatenated list of all unique feature IDs from SEED subsystems
  and literature.  Each line has one field that is the feature ID in KBase format.
* **CDM_ALL_FID_ROLE**: Mapping of all feature IDs to functional roles.  Each line
  has two fields: (1) Feature ID in KBase format (2) List of functional roles
  delimited by separator configuration variable.
* **CDM_OTU_FID_ROLE**: Mapping of feature IDs to functional roles where there is
  one protein from each OTU for each functional role.  Each line has two fields:
  (1) Feature ID in KBase format (2) List of functional roles delimited by
  separator configuration variable.
* **CDM_PROTEIN_FASTA**: Fasta file containing the amino acid sequences for proteins
  in list of OTU feature IDs.
* **CDM_COMPLEX_ROLE**: Mapping of complex IDs to functional roles.  Each line has
  two fields: (1) Complex ID in KBase format (2) List of functional roles for the
  complex delimited by separator configuration variable.
* **CDM_REACTION_COMPLEX**: Mapping of reaction IDs to complex IDs. Each line has two
  fields: (1) Reaction ID in KBase format (2) List of complex IDs in KBase format
  delimited by separator configuration variable.

An additional source is KEGG and includes these intermediate files:

* **KEGG_OTU_FID_ROLE**: Mapping of feature IDs to functional roles.  Each line has
  two fields: (1) Feature ID in KEGG format (2) List of functional roles delimited by
  separator configuration variable.
* **KEGG_PROTEIN_FASTA**: Fasta file containing the amino acid sequences for proteins
  in list of OTU feature IDs.
* **KEGG_COMPLEX_ROLE**: Mapping of complex IDs to functional roles.  Each line has
  two fields: (1) Complex ID in KBase format (2) List of functional roles for the
  complex delimited by separator configuration variable.
* **KEGG_REACTION_COMPLEX**: Mapping of reaction IDs to complex IDs. Each line has two
  fields: (1) Reaction ID in KBase format (2) List of complex IDs in KBase format
  delimited by separator configuration variable.

For building ProbAnno and RxnProbs typed objects, the intermediate files from the
configured sources are merged into these final files:

* **OTU_FID_ROLE**: Mapping of feature IDs to functional roles.  Each line has two
  fields: (1) Feature ID in source format (2) List of functional roles delimited by
  separator configuration variable.
* **PROTEIN_FASTA**: Fasta file containing the amino acid sequences for proteins
  in list of OTU feature IDs.
* **COMPLEX_ROLE**: Mapping of complex IDs to functional roles.  Each line has
  two fields: (1) Complex ID in KBase format (2) List of functional roles for the
  complex delimited by separator configuration variable.
* **REACTION_COMPLEX**: Mapping of reaction IDs to complex IDs. Each line has two
  fields: (1) Reaction ID in KBase format (2) List of complex IDs in KBase format
  delimited by separator configuration variable.

Building the Static Database Files
----------------------------------

The static database files are built in multiple steps from one or more sources.
The Probabilistic Annotation methods use a merged set of files for calculating
likelihoods.  The files are managed with the following commands:

* **pa-gendata-cdm**: Generates intermediate files from KBase Central Data Model.
* **pa-gendata-kegg**: Generates intermediate files from Kyoto Encylopedia of
  Genes and Genomes.
* **pa-gendata**: Merges intermediate files into final files and builds search
  database for configured search program.
* **pa-savedata**: Stores final files in Shock.
* **pa-loaddata**: Retrieves final files from Shock.

Build the static database files by following these steps:

1. Run the pa-gendata-cdm command to generate KBase Central Data Model
   intermediate files.
2. Optionally, run the pa-gendata-kegg command to generate KEGG intermediate
   files (requires local copies of KEGG reaction, enzyme, and organism databases).
3. Run the pa-gendata command to generate final files and search database.


Loading the Static Database Files
---------------------------------

The static database files can be loaded for the server in these ways:

1. Automatically from Shock.  When the load\_data_option
   configuration variable is set to "shock", the server checks Shock to
   see if a new version of each file is available and downloads any new
   or missing files. 
2. When the load\_data_option configuration variable is set to "preload",
   use the files already loaded in the data folder.
