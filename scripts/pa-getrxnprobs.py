import argparse
import traceback
import sys
from biokbase.probabilistic_annotation.Helpers import get_url
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation
from biokbase.workspace.ScriptHelpers import user_workspace

desc1 = '''
NAME
      pa-getrxnprobs -- get a table of reaction probabilities from a RxnProbs object

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Get a table of reaction probabilities from a RxnProbs object.  Each
      reaction in the RxnProbs object is given a single row in the output table.
      Reactions with no complexes linked to them will have no rows in the table.

      The --workspace optional argument specifies the workspace where the
      RxnProbs object is stored.  By default, the user's current workspace as
      set by the ws-workspace command is used to find the object.

      The --version optional argument specifies the version number of the
      RxnProbs object.  By default, the latest version is used.

      The --sort optional argument specifies the field to use as a key for
      sorting the output table.  Valid values are "rxnid" or "probability".  The
      default is to sort using the reaction ID as the key.

      The --url optional argument specifies an alternate URL for the service
      endpoint.

      The --show-error optional argument shows additional detailed information
      when an exception occurs.
'''

desc3 = '''
EXAMPLES
      > pa-getrxnprobs 'kb|g.0.rxnprobs'
      reaction_id   probability   complex_diagnostic   complex_details   putative_GPR

SEE ALSO
      pa-calculate
      ws-workspace

AUTHORS
      Matt Benedict, Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-getrxnprobs', epilog=desc3)
    parser.add_argument('rxnprobs', help='ID of RxnProbs object', action='store', default=None)
    parser.add_argument('-w', '--workspace', help='workspace where RxnProbs object is saved', action='store', dest='rxnprobsws', default=None)
    parser.add_argument('-v', '--version', help='version number of RxnProbs object', action='store', dest='rxnprobsver', type=int, default=None)
    parser.add_argument('--sort', help='field to use as key for sorting output table', action='store', dest='sortField', default='rxnid')
    parser.add_argument('-u', '--url', help='url for service', action='store', dest='url', default=None)
    parser.add_argument('-e', '--show-error', help='show detailed information for an exception', action='store_true', dest='showError', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create input parameters for annotate() function.
    input = dict()
    input['rxnprobs'] = args.rxnprobs
    if args.rxnprobsws is None:
        input['rxnprobs_workspace'] = user_workspace()
    else:
        input['rxnprobs_workspace'] = args.rxnprobsws
    input['rxnprobs_version'] = args.rxnprobsver
    input['sort_field'] = args.sortField
                
    # Create a probabilistic annotation client.
    if args.url is None:
        args.url = get_url()
    paClient = ProbabilisticAnnotation(url=args.url)

    # Get the raw data from the server.
    try:
        output = paClient.get_rxnprobs(input)
    except Exception as e:
        print 'Error getting RxnProbs object: %s' %(e.message)
        if args.showError:
            traceback.print_exc(file=sys.stdout)
        exit(1)

    # Format the data as a table of tab delimited fields.
    print '\t'.join(['reaction_id', 'probability', 'complex_diagnostic', 'complex_details', 'putative_GPR'])
    for rxnprob in output:
        print '%s\t%f\t%s\t%s\t%s' %(rxnprob[0], rxnprob[1], rxnprob[2], rxnprob[3], rxnprob[4]) 

    exit(0)

    
