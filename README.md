Probabilistic Annotation
========================

The purpose of the Probabilistic Annotation service is to provide users with
alternative annotations for genes, each attached to a likelihood score, and to
translate these likelihood scores into likelihood scores for the existence of
reactions in metabolic models.  With the Probabilistic Annotation service:

* Users can quickly assess the quality of an annotation.

* Reaction likelihood computations allow users to estimate the quality of
  metabolic networks generated using the automated reconstruction tools in
  other services.

* Combining reaction likelihoods with gap filling both directly incorporates
  available genetic evidence into the gap filling process and provides putative
  gene annotations automatically, reducing the effort needed to search for
  evidence for gap filled reactions.

Details on the algorithm available in [Likelihood-Based Gene Annotations
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
* **work_folder_path**: Path to work folder containing sub-folders for running service.
  requests.  The intermediate files created by the annotate() and calculate() methods
  are stored in this location.
* **data_folder_path**: Path to data folder containing static database files.  A server
  can use different static database files by using different locations.
* **load_data_option**: Control how static database files are handled when starting
  service. Valid values are "shock" to load static files from Shock or 'preload'
  to use preloaded static database files.
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

