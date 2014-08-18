import argparse
from biokbase.probabilistic_annotation.Helpers import get_url, set_url
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation, ServerError as ProbAnnoServerError

desc1 = '''
NAME
      pa-url -- update or view url of the probabilistic annotation service endpoint

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Display or set the URL endpoint for the probabilistic annotation service.
      If run with no arguments or options, then the current URL is displayed.
      If run with a single argument, the current URL is set to the specified
      URL.  If the specified URL is named default, then the URL is reset to
      the default production URL.  When the --no-check optional argument is
      specified, the validity of the endpoint is not checked.
'''

desc3 = '''
EXAMPLES
      Display the current URL:
      > pa-url
      Current URL: https://kbase.us/services/probabilistic_annotation/

      Use a new URL:
      > pa-url http://localhost:7073
      http://localhost:7073 is valid and running probabilistic_annotation v1.1.0

      Reset to the default URL:
      > pa-url default
      https://kbase.us/services/probabilistic_annotation/ is valid and running probabilistic_annotation v1.1.0

AUTHORS
      Matt Benedict, Mike Mundy
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa_url', epilog=desc3)
    parser.add_argument('newurl', nargs='?', default=None, help='New URL endpoint')
    parser.add_argument('--no-check', help='do not check for a valid URL', action='store_true', dest='noCheck', default=False)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    if args.newurl == None:
        print "Current URL: " + get_url()
    else:
        newurl = set_url(args.newurl)
        if args.noCheck:
            print 'New URL set to: '+newurl
        else:
            try:
                paClient = ProbabilisticAnnotation(url=newurl)
                serverInfo = paClient.version()
                if serverInfo[0] != 'probabilistic_annotation':
                    print newurl+' is not a probabilistic annotation server'
                    exit(1)
                print newurl+' is valid and running %s v%s' %(serverInfo[0], serverInfo[1])
            except ProbAnnoServerError as e:
                print 'Endpoint at %s returned error: %s' %(newurl, e)
                exit(1)
            except Exception as e:
                print 'Could not get a valid response from endpoint at %s: %s' %(newurl, e)
                exit(1)
    exit(0)
    