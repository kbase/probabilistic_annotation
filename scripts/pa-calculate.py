import argparse
import traceback
import sys
from biokbase.probabilistic_annotation.Helpers import get_url
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation
from biokbase.workspace.ScriptHelpers import user_workspace

desc1 = '''
NAME
      pa-calculate -- calculate reaction likelihoods from a probabilistic annotation

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Calculate reaction likelihoods from a probabilistic annotation generated
      by the pa-annotate command.
      
      The results are saved in a RxnProbs typed object and contain putative
      gene annotations (based on a cutoff from the gene most likely to fulfill
      each role associated with the reaction) and likelihood scores.
      
      The RxnProbs object can be used as input to gapfilling a metabolic model
      using the --probrxn option for the fba-gapfill command.  However, if 
      you do this you must make sure that the same template model is used for
      gapfilling and for computing probabilities.  If you want to avoid this
      issue, we recommend using the ProbAnno object instead.
      
      (default is to use all reactions in the biochemistry)
'''

desc3 = '''
EXAMPLES
      Calculate reaction likelihoods from probabilistic annotation of E. coli
      K12 genome:
      > pa-calculate kb|g.0.probanno kb|g.0.rxnprobs

SEE ALSO
      pa-annotate
      pa-getrxnprobs
      fba-gapfill

AUTHORS
      Matt Benedict, Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-calculate', epilog=desc3)
    parser.add_argument('probanno', help='ID of ProbAnno object', action='store', default=None)
    parser.add_argument('rxnprobs', help='ID of RxnProbs object', action='store', default=None)
    parser.add_argument('-w', '--rxnprobsws', help='workspace where RxnProbs object is saved', action='store', dest='rxnprobsws', default=None)
    parser.add_argument('--probannows', help='workspace where ProbAnno object is stored', action='store', dest='probannows', default=None)
    parser.add_argument('-t', '--template', help='ID of ModelTemplate object', action='store', dest='template', default=None)
    parser.add_argument('--templatews', help='workspace where ModelTemplate object is stored', action='store', dest='templatews', default=None)
    parser.add_argument('--url', help='url for service', action='store', dest='url', default=None)
    parser.add_argument('-e', '--show-error', help='show detailed information for an exception', action='store_true', dest='showError', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    # Create input parameters for annotate() function.
    input = dict()
    input['probanno'] = args.probanno
    input['rxnprobs'] = args.rxnprobs
    if args.probannows is None:
        input['probanno_workspace'] = user_workspace()
    else:
        input['probanno_workpace'] = args.probannows
    if args.rxnprobsws is None:
        input['rxnprobs_workspace'] = user_workspace()
    else:
        input['rxnprobs_workspace'] = args.rxnprobsws
    input['template_model'] = args.template
    input['template_workspace'] = args.templatews
                
    # Create a probabilistic annotation client.
    if args.url is None:
        args.url = get_url()
    paClient = ProbabilisticAnnotation(url=args.url)

    # Submit a job to annotate the specified genome.
    try:
        jobid = paClient.calculate(input)
        # python version of print object info
    except Exception as e:
        print 'Error calculating reaction probabilities: %s' %(e.message)
        if args.showError:
            traceback.print_exc(file=sys.stdout)
        exit(1)

    exit(0)
