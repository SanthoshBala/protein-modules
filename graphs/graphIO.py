#! /usr/bin/python

# graphIO.py
# Author: Santhosh Balasubramanian
# Created: March 12, 2013
# Last Modified: March 24, 2013


# Library Imports
from graph_tool.all import *

# Global Imports
from settings import *

# Utility Imports
from common.strings import *


# - - - - - - - - - - GRAPH INPUT - - - - - - - - - - #


# readGraph: Returns graph from <filename>. If <canonical> is true,
# then edge weights represent confidences.
def readGraph(filename, canonical = False, weighted = True, directed = False):
    # Open File
    inFile = open(filename, 'r')

    # Create graph-tool Graph
    if directed:
        graph = Graph(directed=True)
    else:
        graph = Graph(directed=False)

    # Add vertex property for gene_id
    graph.vertex_properties['gene_id'] = graph.new_vertex_property('string')
    graph.edge_properties['confidence'] = graph.new_edge_property('double')

    # Optionally skip a line
    skipLine = True

    # Note that we must keep track of which genes already have a 
    # corresponding node in the graph. We could just use the built in
    # PropertyMap of the graph, but this is just an array of length N,
    # where array[i] is the value for the ith vertex. Lookup is therefore
    # extremely slow (linear), so we should use a python dict (constant)
    geneDict = dict()
    geneIDProperty = graph.vertex_properties['gene_id']

    # Read all lines in network file
    for line in inFile:
        # Ignore header line
        if skipLine:
            skipLine = False
            continue
        
        # Get Fields from Genemania File
        lineFields = parseTabSeparatedLine(line)
        geneAName = lineFields[0]
        geneBName = lineFields[1]
        confStr = lineFields[2]
        if canonical:
            confFloat = float(confStr)
        else:
            confFloat = 1.0 - float(confStr)
        
        if confFloat == 0.0:
            continue
        
        # Get geneA's vertex index if it exists
        if geneAName in geneDict:
            geneANode = geneDict[geneAName]
        else:
            geneANode = graph.add_vertex()
            geneIDProperty[geneANode] = geneAName
            geneDict.update( { geneAName : geneANode } )
             
        # Get geneB's vertex if it exists
        if geneBName in geneDict:
            geneBNode = geneDict[geneBName]
        else:
            geneBNode = graph.add_vertex()
            geneIDProperty[geneBNode] = geneBName
            geneDict.update( { geneBName : geneBNode } )

        # Create edge between geneANode and geneBNode
        newEdge = graph.add_edge(geneANode, geneBNode)
        if weighted:
            graph.edge_properties['confidence'][newEdge] = confFloat
        else:
            graph.edge_properties['confidence'][newEdge] = 1.0


    # Close file
    inFile.close()

    # Return graph object
    return graph

# getTissueSubgraph: Reads PPI subgraph from file.
def getTissueSubgraph(tissue):
    interactionType = DESIRED_INTERACTION_TYPE

    # Construct File Path
    filename = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, interactionType)
    filepath = PATH_TO_TISSUE_SUBGRAPHS + filename

    # Read Graph to Object
    graph = readGraph(filepath)
    return graph


# getUnweightedTissueSubgraph: Reads PPI subgraph from file.
def getUnweightedTissueSubgraph(tissue):
    interactionType = DESIRED_INTERACTION_TYPE

    # Construct File Path
    filename = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, interactionType)
    filepath = PATH_TO_TISSUE_SUBGRAPHS + filename

    # Read Graph to Object
    graph = readGraph(filepath, weighted = False)
    return graph


# - - - - - - - - - - GRAPH OUTPUT - - - - - - - - - - #


# writeGraph: Writes <graph> to <outFilePath>.
#  If <canonical> is true, then edge weights represent confidences.
def writeGraph(graph, outFilePath, canonical = False, weighted = True):
    # Open Output File
    outFile = open(outFilePath, 'w')

    # Get Gene Names and Interaction Confidences
    geneIDProp = graph.vertex_properties['gene_id']
    confidenceProp = graph.edge_properties['confidence']

    outFile.write('Gene_A\tGene_B\tWeight\n')

    # Iterate Over All Edges in Graph
    for edge in graph.edges():
        source = edge.source()
        dest = edge.target()

        # Get Line Fields
        sourceName = geneIDProp[source]
        destName = geneIDProp[dest]
        if canonical:
            confStr = str(confidenceProp[edge])
        else:
            confStr = str(1.0 - confidenceProp[edge])

        if not weighted:
            confStr = '1.0'

        # Write to File
        outFile.write('%s\t%s\t%s\n' % (sourceName, destName, confStr))

    # Close File
    outFile.close()


# - - - - - - - - - - FORMAT CONVERSION - - - - - - - - - - #


# canonicalizeGenemaniaGraph: Converts Genemania format into 
# format recognizable by SPICi, with confidence in 3rd column.
def canonicalizeGenemaniaGraph(filename):
    # Get Output Path
    outFilename = 'canon.%s' % filename
    outFilepath = PATH_TO_CANONICAL_TISSUE_SUBGRAPHS + outFilename

    # Read <filename> into Graph
    graph = readGraph(filename)

    # Write to Disk
    writeGraph(graph, outFilepath, canonical = True)

    return

# convertCanonicalGraphToSIF: Converts canonical format to 
# SIF format, for use with cytoscape 3.0
def convertCanonicalGraphToSIF(inFilepath, inFilename):
    # Open Input File
    inFile = open(inFilepath + inFilename, 'r')

    # Open Output File
    outFilename = inFilename + '.sif'
    outFile = open(inFilepath + outFilename, 'w')

    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        outFile.write("%s\t%s\t%s\n" % (lineFields[0], lineFields[2], lineFields[1]))

    # Close Files
    inFile.close()
    outFile.close()

    return
