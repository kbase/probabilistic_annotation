import argparse
import traceback
import sys
from biokbase.probabilistic_annotation.Helpers import get_url
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation
from biokbase.workspace.ScriptHelpers import user_workspace

desc1 = '''
NAME
      pa-getprobanno -- get a table of probabilistic annotations

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Get a table of annotation probabilities from a ProbAnno object.  Each
      gene-annotation pair is given its own row in the table. An annotation is
      a set of roles delimited by the separator '///'.

      The --workspace optional argument specifies the workspace where the
      ProbAnno object is stored.  By default, the user's current workspace as
      set by the ws-workspace command is used to find the object.

      The --roles optional argument gets gene-role pairs instead where the
      probability of the role is computed as the sum of the probabilities of
      annotations containing it.

      The --url optional argument specifies an alternate URL for the service
      endpoint.

      The --show-error optional argument shows additional detailed information
      when an exception occurs.
'''

desc3 = '''
EXAMPLES
      > pa-getprobanno 'kb|g.0.probanno'
      gene    annotation   likelihood

      > pa-getprobanno --roles 'kb|g.0.probanno'
      gene    role    likelihood

SEE ALSO
      pa-annotate
      ws-workspace

AUTHORS
      Matt Benedict, Mike Mundy 
'''

if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa-getprobanno', epilog=desc3)
    parser.add_argument('probanno', help='ID of ProbAnno object', action='store', default=None)
    parser.add_argument('-w', '--workspace', help='workspace where ProbAnno object is saved', action='store', dest='probannows', default=None)
    parser.add_argument('-r', '--roles', help='Print role likelihoods instead of annotation likelihoods', action='store_true', dest='roles', default=False)
    parser.add_argument('-u', '--url', help='url for service', action='store', dest='url', default=None)
    parser.add_argument('-e', '--show-error', help='show detailed information for an exception', action='store_true', dest='showError', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()

    # Create input parameters for annotate() function.
    input = dict()
    input['probanno'] = args.probanno
    if args.probannows is None:
        input['probanno_workspace'] = user_workspace()
    else:
        input['probanno_workspace'] = args.probannows
                
    # Create a probabilistic annotation client.
    if args.url is None:
        args.url = get_url()
    paClient = ProbabilisticAnnotation(url=args.url)

    # Get the raw data from the server.
    try:
        output = paClient.get_probanno(input)
    except Exception as e:
        print 'Error getting ProbAnno object: %s' %(e.message)
        if args.showError:
            traceback.print_exc(file=sys.stdout)
        exit(1)

    # Format the data as a table of tab delimited fields.
    if args.roles:
        print '\t'.join(['gene', 'role', 'likelihood'])
    else:
        print '\t'.join(['gene', 'annotation', 'likelihood'])
    for gene in output:
        roleToProb = dict()
        for roleprob in output[gene]:
            if args.roles:
                # Extract the roles from the roleset string and add them up if there are duplicates.
                roles = roleprob[0].split('///')
                for role in roles:
                    if role in roleToProb:
                        roleToProb[role] += roleprob[1]
                    else:
                        roleToProb[role] = roleprob[1]
            else:
                roleToProb[roleprob[0]] = roleprob[1]
        for role in roleToProb:
            print '%s\t%s\t%f' %(gene, role, roleToProb[role])

    exit(0)

    
