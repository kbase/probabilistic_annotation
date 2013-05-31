import sys
import unittest
import subprocess
import time
import os


# Used for cleanup after the test
from biokbase.workspaceService.Client import *

class TestAnnotateScript(unittest.TestCase):    
    def setUp(self):
        self.cmd = os.environ["KB_TOP"] + "/bin/pa-annotate"
        self.url = "http://localhost:7073"
        self.wsurl = "http://kbase.us/services/workspace/"

        self.genomeid = "TEST_SMALL_GENOME"
        self.genomews = "ProbAnnoTest"
        self.probannoid = "TEST_SMALL_GENOME.probanno"
        self.probannows = "ProbAnnoTest"

        # When the auth service is implemented in Python I'll replace this.
        self.tokenfile = os.path.join(os.path.expanduser("~"), ".kbase_auth")
        self.auth = "".join( [ line.strip("\r\n") for line in open(self.tokenfile) ] )

        # This must be at leas long enough for the command to run on the genes in the test object.
        self.runtime = 30

        # FIXME: This option doesn't quite work... it still creates a new version rather than just removing
        # the old one and writing a new one. Either that or that's not the intention.
        self.overwrite = True

    '''
    Test inputs and option processing
    '''
    def test_help(self):
        '''Run pa-annotate --help and verify that the major sections in the help text are present'''
        
        args = [ self.cmd, "--help" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find("NAME"), -1)
        self.assertNotEqual(so.find("SYNOPSIS"), -1)
        self.assertNotEqual(so.find("DESCRIPTION"), -1)
        self.assertNotEqual(so.find("EXAMPLES"), -1)
        self.assertEqual(se, '')        
    def test_badOption(self):
        '''Run pa-annotate with a bad option and verify that the error message is returned.'''
        
        args = [ self.cmd, "kb|g.8622", "kb|g.8622", "--chia" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 255)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Unknown option:"), -1)
    def test_missingOptionValue(self):
        '''Run pa-annotate with a missing option value and verify that the error message is returned.'''
        
        args = [ self.cmd, "kb|g.8622", "kb|g.8622", "--genomews" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 255)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Option genomews requires an argument"), -1)
    def test_missingArg(self):
        '''Run pa-annotate with a missing argument and verify that the error message is returned.'''
        
        args = [ self.cmd, "kb|g.8622", "--genomews", self.genomews ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Required arguments are missing"), -1)
    def test_run(self):
        ''' Run pa-annotate on a valid genome object and verify that the job runs and returns a result in the expected time'''

        args = [ self.cmd, self.genomeid, self.probannoid, "--genomews", self.genomews, "-w", self.probannows ]
        if self.overwrite:
            args.append("-o")

        wsClient = workspaceService(self.wsurl)

        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (sc, se) = proc.communicate()

        # Running the command must not fail (e.g. due to missing arguments).
        self.assertEqual(proc.returncode, 0)

        # Allow time for the command to run
        time.sleep(30)

        # Lets grab the object
        wsClientArgs = { "workspace": self.probannows,
                         "type": "ProbAnno",
                         "id": self.probannoid,
                         "auth" : self.auth,
                         }

        # Did we successfully make a probanno object?
        try:
            probanno_object = wsClient.get_object(wsClientArgs)
        except ServerError:
            self.fail(msg="The expected object %s did not get created in the output workspace!\n" %(outputobj))
        except:
            self.fail(str(sys.exc_info()[0]))

if __name__ == '__main__':
    # Since I can't get this from python I think...
    unittest.main()
