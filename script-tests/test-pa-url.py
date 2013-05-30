import sys
import unittest
import subprocess
from os import environ

class TestUrlScript(unittest.TestCase):
    
    def setUp(self):
        self.cmd = environ["KB_TOP"] + "/bin/pa-url"
        self.url = "http://localhost:7073"
        
    def test_current(self):
        '''Run pa-url and verify that the current url is returned.'''
        args = [ self.cmd ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find("http://"), -1)
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
        
        args = [ self.cmd, "default" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        self.assertNotEqual(so.find("http://kbase.us/services/probabilistic_annotation/"), -1)
        self.assertEqual(se, '')
        
    def test_seturl(self):
        '''Run pa-url newurl and verify that the new URL is returned.'''
        
        args = [ self.cmd, self.url ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 0)
        lines = so.split("\n")
        self.assertEqual(lines[1], self.url)
        self.assertEqual(se, '')
        
    def test_badarg(self):
        '''Run pa-url with a bad argument and verify that the error message is returned.'''
        
        args = [ self.cmd, "--chia" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 255)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Unknown option:"), -1)

if __name__ == '__main__':
    unittest.main()
