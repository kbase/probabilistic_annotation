import sys
import unittest
import subprocess
import os
from biokbase.probabilistic_annotation.Helpers import get_config

class TestUrlScript(unittest.TestCase):
    
    def setUp(self):
        self.cmd = os.path.join(os.environ["KB_TOP"], "bin/pa-url")
        self.config = get_config(os.environ["KB_TEST_CONFIG"])
        
    def test_current(self):
        '''Run pa-url and verify that the current url is returned.'''
        args = [ self.cmd ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find("Current URL:"), -1)
        self.assertEqual(se, '')
        lines = so.split("\n")
        self.url = lines[1]
       
    def test_help(self):
        '''Run pa-url --help and verify that the major sections in the help text are present'''
        
        args = [ self.cmd, "--help" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find("NAME"), -1)
        self.assertNotEqual(so.find("SYNOPSIS"), -1)
        self.assertNotEqual(so.find("DESCRIPTION"), -1)
        self.assertNotEqual(so.find("EXAMPLES"), -1)
        self.assertEqual(se, '')
        
    def test_default(self):
        '''Run pa-url default and verify that the default URL is returned.'''
        
        args = [ self.cmd, "default", "--no-check" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find("https://kbase.us/services/probabilistic_annotation/"), -1)
        self.assertEqual(se, '')
        
    def test_seturl(self):
        '''Run pa-url newurl and verify that the new URL is returned.'''
        
        args = [ self.cmd, self.config["probanno_url"] ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find(self.config["probanno_url"]), -1)
        self.assertEqual(se, '')
        
    def test_badarg(self):
        '''Run pa-url with a bad argument and verify that the error message is returned.'''
        
        args = [ self.cmd, "--chia" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertNotEqual(proc.returncode, 0)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("unrecognized arguments:"), -1)

if __name__ == '__main__':
    unittest.main()
