#! /usr/bin/python

# shell.py
# Author: Santhosh Balasubramanian
# Created: March 24, 2013
# Last Modified: March 24, 2013


# Python Imports
import sys
import subprocess

# Global Imports
from settings import *


# mountCodeDir: Add '/code/' to PYTHONPATH.
def mountCodeDir():
    sys.path.append(PATH_TO_CODE)

# getFileLineCount: Get # of lines in <filepath>.
def getFileLineCount(filepath):
    # Construct 'wc -l' Command
    command = ['wc', '-l', filepath]

    # Get Line Count
    output = subprocess.check_output(command)
    lineCount = int(output.split()[0])

    return lineCount

