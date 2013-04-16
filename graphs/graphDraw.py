#! /usr/bin/python

# graphDraw.py
# Author: Santhosh Balasubramanian
# Created: April 15, 2013
# Last Modified: April 15, 2013


# Library Imports
from graph_tool.all import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphUtil import *

# Module Imports
from modules.moduleUtil import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

PATH_TO_REPORT_IMAGES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/'
MODULE_DRAWING_BASE_FILENAME = '%s-graph.pdf'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# drawModule: Draws <moduleID> graph to pdf.
def drawModule(moduleID):

    # Get Module Graph
    moduleGraph = getModuleGraph(moduleID)
    
    outFileName = MODULE_DRAWING_BASE_FILENAME % moduleID
    outFile = PATH_TO_REPORT_IMAGES + outFileName

    # Get Vertex Positions

    # Get Vertex Text
    geneIDProp = moduleGraph.vertex_properties['gene_id']

    # Draw Module
    graph_tool.draw.graph_draw(moduleGraph, 
                               output = outFile, vertex_text = geneIDProp)

    return

# drawGraph
def drawGraph():
    g = graph_tool.generation.price_network(10)
    d = g.degree_property_map('out')

    outFileName = ''
    graph_tool.draw.graph_draw(g, output='/home/santhosh/graph.draw.test.pdf')
    graph_tool.draw.graphviz_draw(g,output='/home/santhosh/graphviz.draw.test.pdf')
