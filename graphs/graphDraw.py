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
from graphs.graphIO import *
from graphs.graphUtil import *

# Module Imports
from modules.moduleUtil import *
from modules.moduleTopology import *

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
                                      ecolor = '#000000',
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
                                      ecolor = '#000000',
                                      vprops = { 'label' : geneIDProp,
                                                 'fontcolor' : '#ffffff',
                                                 'fontname' : 'Palatino',
                                                 'penwidth' : 1.0}, 
                                      output = outFile)

    return


# drawTissueModuleTopology
def drawTissueModuleTopology(tissue):
    # Get Tissue Subgraph
    graph = getTissueSubgraph(tissue)

    # Get Module Topology
    topology = getModuleTopology(graph, initialize = False)

    # Get Output File
    outFilePath = '/home/santhosh/Dropbox/Thesis/report/final_report/images/'
    outFileName = '%s-module-topology.pdf' % tissue
    outFile = outFilePath + outFileName

    # Draw Topology
    heightProp = topology.vertex_properties['height']
    graph_tool.draw.graphviz_draw(topology, 
                                  overlap = False, 
                                  vcolor = heightProp,
                                  output = outFile
                                  )
    

    return


# drawRavaszBarabasiModel
def drawRavaszBarabasiModel():
    # Specify Vertex Colors
    colorDict = {}
    for i in range(1, 65):
        if i < 5:
            colorDict.update( { str(i) : '#3333FF'} )
        elif i < 17:
            colorDict.update( { str(i) : '#007C00'} )
        else:
            colorDict.update( { str(i) : '#CC0000' } )

    # Read Round One Graph
    INPUT_BASE_FILENAME = 'ravasz.barabasi.%d'
    OUTPUT_BASE_FILENAME = 'ravasz-barabasi-%d.pdf'
    inFilePath = '/home/santhosh/workspace/thesis/data/analysis/module_topology/'
    outFilePath = '/home/santhosh/Dropbox/Thesis/report/final_report/images/'
    
    for i in range(1,3):
        # Get Input File Name
        inFileName = INPUT_BASE_FILENAME % i
        inFile = inFilePath + inFileName
        
        # Get Output File Name
        outFileName = OUTPUT_BASE_FILENAME % i
        outFile = outFilePath + outFileName
        
        graph = readGraph(inFile, canonical = True)
        # Color Vertices
        graph.vertex_properties['color'] = graph.new_vertex_property('string')
        colorProp = graph.vertex_properties['color']
        geneIDProp = graph.vertex_properties['gene_id']
        for vertex in graph.vertices():
            geneID = geneIDProp[vertex]
            color = colorDict.get(geneID)
            colorProp[vertex] = color

        graph_tool.draw.graphviz_draw(graph, 
                                      overlap = False, 
                                      vcolor = colorProp,
                                      layout = 'twopi',
                                      ecolor = '#000000',
                                      vprops = {'fixedsize' : True, },
                                      gprops = { 'rankdir' : 'BT'
                                                 },
                                      output = outFile
                                      )

    return


# - - - - - - - - - - GRAPH MODEL - - - - - - - - - - #

# drawRavaszGraphModuleTopology
def drawRavaszGraphModuleTopology():
    
    colorDict = {
        'A' : '#CC0000',
        'B' : '#CC0000',
        'C' : '#CC0000',
        'D' : '#993399',
        'E' : '#007C00',
        'F' : '#007C00',
        'G' : '#007C00',
        'H' : '#3333FF',
        'I' : '#3333FF',
        'J' : '#3333FF',
        'K' : '#3333FF',
        'mod:A' : '#CC0000',
        'mod:E' : '#007C00',
        'mod:H' : '#3333FF',
        'mod:D' : '#999999',
        'mod:mod:A' : '#999999',
        }

    # ROUND 1

    # Read Graph
    inFilePath = '/home/santhosh/workspace/thesis/data/analysis/module_topology/'
    inFileName = 'ravasz.sample.graph'
    graph = readGraph(inFilePath + inFileName, canonical = True)

    # Add Colors to Graph
    graph.vertex_properties['color'] = graph.new_vertex_property('string')
    colorProp = graph.vertex_properties['color']
    geneIDProp = graph.vertex_properties['gene_id']
    for vertex in graph.vertices():
        geneID = geneIDProp[vertex]
        color = colorDict.get(geneID)
        colorProp[vertex] = color
    
    # Get Output File Name
    outFilePath = PATH_TO_REPORT_IMAGES
    outFileName = 'ravasz-graph.pdf'
    outFile = outFilePath + outFileName

    # Draw Graph
    penWidth = 2.0
    graph_tool.draw.graphviz_draw(graph, 
                                  overlap = False, 
                                  vsize = 0.5,
                                  vcolor = colorProp,
                                  ecolor = '#000000',
                                  vprops = { 'label' : geneIDProp,
                                             'fontsize' : 12,
                                             'fontcolor' : '#ffffff',
                                             'fontname' : 'Palatino',
                                             'penwidth' : penWidth,
                                             'fixedsize' : True }, 
                                  output = outFile)
    
    # Get Module Topology
    topology = getModuleTopology(graph)
    geneIDProp = topology.vertex_properties['gene_id']
    topology.vertex_properties['color'] = topology.new_vertex_property('string')
    colorProp = topology.vertex_properties['color']

    topologyLabelDict = {
        'A' : 'A',
        'B' : 'B',
        'C' : 'C',
        'D' : 'D',
        'E' : 'E',
        'F' : 'F',
        'G' : 'G',
        'H' : 'H',
        'I' : 'I',
        'J' : 'J',
        'K' : 'K',
        'mod:A' : 'M1',
        'mod:E' : 'M2',
        'mod:H' : 'M3',
        'mod:mod:A' : 'M4',
        'mod:D' : 'M4',
        }

    # Relabel Vertices and Reassign Colors
    for vertex in topology.vertices():
        geneID = geneIDProp[vertex]
        color = colorDict.get(geneID)
        colorProp[vertex] = color

        newGeneID = topologyLabelDict.get(geneID)
        geneIDProp[vertex] = newGeneID

    # Get Output File Name
    outFileName = 'ravasz-module-topology.pdf'
    outFile = outFilePath + outFileName

    graph_tool.draw.graphviz_draw(topology, 
                                  overlap = False, 
                                  vsize = 0.5,
                                  vcolor = colorProp,
                                  layout = 'dot',
                                  ecolor = '#000000',
                                  vprops = { 'label' : geneIDProp,
                                             'fontsize' : 12,
                                             'fontcolor' : '#ffffff',
                                             'fontname' : 'Palatino',
                                             'penwidth' : penWidth,
                                             'fixedsize' : True, },
                                  gprops = { 'rankdir' : 'BT'
                                             },
                                  output = outFile)
    

    return

