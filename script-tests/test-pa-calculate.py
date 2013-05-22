import sys
import unittest
import subprocess
from os import environ

class TestCalculateScript(unittest.TestCase):
    
    def setUp(self):
        self.cmd = environ["TARGET"] + "/bin/pa-calculate"
        self.url = "http://localhost:7073"
        
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

if __name__ == '__main__':
    unittest.main()