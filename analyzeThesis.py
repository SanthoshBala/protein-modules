#! /usr/bin/python

# analyzeThesis.py
# Author: Santhosh Balasubramanian
# Date: January 23, 2013

import sys
from settings import *
from common import *

from expression.generateProbeGeneMap import *
from expression.generateSampleProbeMap import *
from expression.generateTissueProbeMap import *
from expression.generateTissueGeneMap import *
from graphs.generateTissueSubgraph import *

OUTPUT_FILENAME = 'analysis.output'

def analyzeThesis():
    
    # Open Output File
    f = open(PATH_TO_ANALYSIS + OUTPUT_FILENAME, 'w')

    # Analyze Probe:Gene Map
    f.write('Analyzing Probe:Gene Map...\n')
    analyzeProbeGeneMap()
    f.write('Success: Analyzed Probe:Gene Map\n')
    f.flush()

    f.write('Analyzing Probe:Gene Map -- Verbose...\n')
    analyzeProbeGeneMap(verbose=True)
    f.write('Success: Analyzed Probe:Gene Map -- Verbose\n')
    f.flush()

    f.write('Analyzing Extended Probe:Gene Map...\n')
    analyzeProbeGeneMap(extended=True)
    f.write('Success: Analyzed Extended Probe:Gene Map\n')
    f.flush()

    f.write('Analyzing Extended Probe:Gene Map -- Verbose...\n')
    analyzeProbeGeneMap(extended=True, verbose=True)
    f.write('Success: Analyzed Extended Probe:Gene Map -- Verbose\n')
    f.flush()

    # Analyze Sample:Probe Map
    #f.write('Analyzing Sample:Probe Map...\n')
    #f.write('Success: Analyzed Sample:Probe Map\n')
    #f.flush()

    # Analyze Tissue:Probe Map
    f.write('Analyzing Tissue:Probe Map...\n')
    analyzeTissueProbeMap(verbose=True)
    f.write('Success: Analyzed Tissue:Probe Map\n')
    f.flush()

    
    # Analyze Tissue:Gene Map
    f.write('Analyzing Tissue:Gene Map...\n')
    analyzeTissueGeneMap()
    f.write('Success: Analyzed Tissue:Gene Map\n')
    f.flush()

    f.write('Analyzing Tissue:Gene Map -- Verbose...\n')
    analyzeTissueGeneMap(verbose=True)
    f.write('Success: Analyzed Tissue:Gene Map -- Verbose\n')
    f.flush()

    # Analyze Tissue Subgraphs
#    f.write('Analyzing Tissue Subgraphs...\n')
#    analyzeAllSubgraphs()
#    f.write('Success: Analyzed Tissue Subgraphs\n')
#    f.flush()

analyzeThesis()
