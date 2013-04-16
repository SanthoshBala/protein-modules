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
def drawModule(moduleID, symbol = False):

    # Get Module Graph
    moduleGraph = getModuleGraph(moduleID, symbol = symbol)
    
    outFileName = MODULE_DRAWING_BASE_FILENAME % moduleID
    outFile = PATH_TO_REPORT_IMAGES + outFileName

    # Get Vertex Positions

    # Get Vertex Text
    geneIDProp = moduleGraph.vertex_properties['gene_id']

    # Draw Module
    if symbol:
        sepVal = 0.25
        penWidth = 3.0
        graph_tool.draw.graphviz_draw(moduleGraph, 
                                      overlap = False, 
#                                      sep = sepVal,
                                      vsize = 0.75,
                                      vcolor = '#993399',
                                      ecolor = '#999999',
                                      vprops = { 'label' : geneIDProp,
                                                 'fontsize' : 12,
                                                 'fontcolor' : '#ffffff',
                                                 'fontname' : 'Palatino',
                                                 'penwidth' : penWidth,
                                                 'fixedsize' : True }, 
                                      output = outFile)
    else:
        graph_tool.draw.graphviz_draw(moduleGraph, 
                                      overlap = False, 
                                      vcolor = '#993399',
                                      ecolor = '#999999',
                                      vprops = { 'label' : geneIDProp,
                                                 'fontcolor' : '#ffffff',
                                                 'fontname' : 'Palatino',
                                                 'penwidth' : 1.0}, 
                                      output = outFile)

    return

# drawGraph
def drawGraph():
    g = graph_tool.generation.price_network(10)
    d = g.degree_property_map('out')

    outFileName = ''
    graph_tool.draw.graph_draw(g, output='/home/santhosh/graph.draw.test.pdf')
    graph_tool.draw.graphviz_draw(g,output='/home/santhosh/graphviz.draw.test.pdf')
