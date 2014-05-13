import argparse
import traceback
import sys
from biokbase.probabilistic_annotation.Helpers import get_url
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation
from biokbase.workspace.ScriptHelpers import user_workspace

desc1 = '''
NAME
      pa-annotate -- generate probabilistic annotation for a genome

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Generate alternative annotations for every gene in a genome together
      with their likelihoods.  The current method for calculating likelihoods
      is based on similarity (BLAST) to genes in subsystems and genes with
      literature evidence.
      
      This command takes a significant amount of time to run (since it has to
      run BLAST against a large database), so it is placed on a queue and 
      returns a job ID.  Use the pa-checkjob command to see if your job has
      finished.  When it is done the results are saved in a ProbAnno typed
      object with the specified ID.
      
      The ProbAnno object can be used as input to gapfilling a metabolic model
      using the --probanno option for the fba-gapfill command.
'''

desc3 = '''
EXAMPLES
      Generate probabilistic annotation for E. coli K12 genome:
      > pa-annotate kb|g.0.genome kb|g.0.probanno

SEE ALSO
      pa-calculate
      pa-checkjob
      pa-getprobanno
      fba-gapfill

AUTHORS
      Matt Benedict, Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-annotate', epilog=desc3)
    parser.add_argument('genome', help='ID of Genome object', action='store', default=None)
    parser.add_argument('probanno', help='ID of ProbAnno object', action='store', default=None)
    parser.add_argument('-w', '--probannows', help='workspace where ProbAnno object is saved', action='store', dest='probannows', default=None)
    parser.add_argument('--genomews', help='workspace where Genome object is saved', action='store', dest='genomews', default=None)
    parser.add_argument('--url', help='url for service', action='store', dest='url', default=None)
    parser.add_argument('-e', '--show-error', help='show detailed information for an exception', action='store_true', dest='showError', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    # Create input parameters for annotate() function.
    input = dict()
    input['genome'] = args.genome
    input['probanno'] = args.probanno
    if args.genomews is None:
        input['genome_workspace'] = user_workspace()
    else:
        input['genome_workspace'] = args.genomews
    if args.probannows is None:
        input['probanno_workspace'] = user_workspace()
    else:
        input['probanno_workspace'] = args.probannows
                
    # Create a probabilistic annotation client.
    if args.url is None:
        args.url = get_url()
    paClient = ProbabilisticAnnotation(url=args.url)

    # Submit a job to annotate the specified genome.
    try:
        jobid = paClient.annotate(input)
        print 'Probabilistic annotation job '+jobid+' successfully submitted'
    except Exception as e:
        print 'Error starting job: %s' %(e.message)
        if args.showError:
            traceback.print_exc(file=sys.stdout)
        exit(1)

    exit(0)
