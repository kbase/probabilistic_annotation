import argparse
from biokbase.probabilistic_annotation.Helpers import get_url, set_url

desc1 = '''
NAME
      pa-url -- update or view url of the probabilistic annotation service endpoint

SYNOPSIS      
'''

desc2 = '''
DESCRIPTION
      Display or set the URL endpoint for the probabilistic annotation service.
      If run with no arguments or options, then the current URL is displayed.
      If run with a single argument, the current URL will be switched to the
      specified URL.  If the specified URL is named default, then the URL is
      reset to the default production URL.
'''

desc3 = '''
EXAMPLES
      Display the current URL:
      > pa-url
      Current URL: https://kbase.us/services/probabilistic_annotation/

      Use a new URL:
      > pa-url http://localhost:7102
      New URL set to: http://localhost:7102

      Reset to the default URL:
      > pa-url default
      New URL set to: https://kbase.us/services/probabilistic_annotation/

AUTHORS
      Matt Benedict, Mike Mundy
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='pa_url', epilog=desc3)
    parser.add_argument('newurl', nargs='?', default=None, help='New URL endpoint')
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    if args.newurl == None:
        print "Current URL: " + get_url()
    else:
        print "New URL set to: " + set_url(args.newurl)
    exit(0)
    