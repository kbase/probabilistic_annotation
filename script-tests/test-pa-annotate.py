import sys
import unittest
import subprocess
from os import environ

class TestAnnotateScript(unittest.TestCase):
    
    def setUp(self):
        self.cmd = environ["TARGET"] + "/bin/pa-annotate"
        self.url = "http://localhost:7073"
        
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
        
        args = [ self.cmd, "kb|g.8622", "--genomews", "chialab" ]
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (so, se) = proc.communicate()
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(so, '')
        self.assertNotEqual(se.find("Required arguments are missing"), -1)

if __name__ == '__main__':
    unittest.main()
