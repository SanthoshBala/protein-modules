#! /usr/bin/python

# strings.py
# Author: Santhosh Balasubramanian
# Created: March 24, 2013
# Last Modified: March 24, 2013


# parseCommaSeparatedLine: Returns fields from comma-separated <line>.
def parseCommaSeparatedLine(line):
    # Split on Comma
    lineCommaSplit = line.split('","')

    # Strip Extra \t and \n
    lineFields = [s.strip('"\t\n') for s in lineCommaSplit]

    return lineFields

# parseTabSeparatedLine: Returns fields from tab-separated <line>.
def parseTabSeparatedLine(line):
    # Split on Tab
    lineTabSplit = line.split('\t')

    # Strip Extra \t
    lineFields = [s.strip('\n') for s in lineTabSplit]

    return lineFields


