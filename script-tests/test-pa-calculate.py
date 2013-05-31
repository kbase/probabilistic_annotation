import sys
import unittest
import subprocess
import os

from biokbase.workspaceService.Client import *

class TestCalculateScript(unittest.TestCase):    
    def setUp(self):
        # FIXME - we really should be getting a lot of these from the config file.
        self.cmd = os.environ["KB_TOP"] + "/bin/pa-calculate"
        self.wsurl = "http://kbase.us/services/workspace/"

        self.probannoid = "TEST_SMALL_GENOME.probanno"
        self.probannows = "ProbAnnoTest"
        self.rxnprobsid = "TEST_SMALL_GENOME.rxnprobs"
        self.rxnprobsws = "ProbAnnoTest"

        self.overwrite = True

        self.tokenfile = os.path.join(os.path.expanduser("~"), ".kbase_auth")
        self.auth = "".join( [ line.strip("\r\n") for line in open(self.tokenfile) ] )
    def test_help(self):
        '''Run pa-calculate --help and verify that the major sections in the help text are present'''
        
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
        '''Run pa-calculate with a bad option and verify that the error message is returned.'''
        
        args = [ self.cmd, "kb|g.8622", "--chia" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 255)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Unknown option:"), -1)

    def test_missingOptionValue(self):
        '''Run pa-calculate with a missing option value and verify that the error message is returned.'''
        
        args = [ self.cmd, "kb|g.8622", "--probannows" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 255)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Option probannows requires an argument"), -1)

    def test_missingArg(self):
        '''Run pa-calculate with a missing argument and verify that the error message is returned.'''
        
        args = [ self.cmd, "--probannows", "chialab" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Required arguments are missing"), -1)

    def test_run(self):
        ''' Run pa-calculate on the specified probanno object and verify that a RxnProbs object is created. '''

        args = [ self.cmd, self.probannoid, self.rxnprobsid, "--probannows", self.probannows, "-w", self.rxnprobsws ]
        # Apparently we haven't implemented this yet
#        if self.overwrite:
#            args.append("-o")

        wsClient = workspaceService(self.wsurl)

        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (sc, se) = proc.communicate()

        # Running the command must not fail (e.g. due to missing arguments).
        self.assertEqual(proc.returncode, 0)

        # Note - we don't need to sleep this time since we aren't sending this to the job queue.
        # If we decide to do so we'll need to sleep here.
        # Since we're not it will be done by the time it gets here so Lets grab the object from the
        # workspace.
        wsClientArgs = { "workspace": self.probannows,
                         "type": "RxnProbs",
                         "id": self.rxnprobsid,
                         "auth" : self.auth,
                         }

        # Did we successfully make a rxnprobs object?
        try:
            probanno_object = wsClient.get_object(wsClientArgs)
        except ServerError:
            self.fail(msg="The expected object %s did not get created in the output workspace!\n" %(outputobj))
        except:
            self.fail(str(sys.exc_info()[0]))


if __name__ == '__main__':
    unittest.main()
