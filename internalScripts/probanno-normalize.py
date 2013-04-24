#!/usr/bin/python

import optparse, sys
from biokbase.probabilistic_annotation.MyImpl import normalize

usage="%prog [options]"
description="""Normalize metabolite weights"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--model", help="ID of model object", action="store", dest="model", default=None)
parser.add_option("-q", "--modelws", help="ID of genome workspace", action="store", dest="model_workspace", default=None)
parser.add_option("-d", "--debug", help="Save generated data for debug purposes", action="store", dest="debug", default=False)
parser.add_option("-v", "--verbose", help="Print verbose output", action="store", dest="verbose", default=False)
parser.add_option("-a", "--auth", help="Auth token", action="store", dest="auth", default=None)
parser.add_option("-b", "--absval", help="Absolute value", action="store", dest="absval", default=False)
(options, args) = parser.parse_args()

success = normalize(options)
if success:
    exit(0)
exit(1)
