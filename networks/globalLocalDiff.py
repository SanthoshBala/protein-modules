#! /usr/bin/python

# Santhosh Balasubramanian
# February 27, 2013

from settings import *
from graphUtil import *
from moduleUtil import *

IGNORE_INHERITED_DIFF = False

def compareGlobalAndLocalGraphs():
    # Get Global Graph
    globalGraph = readGenemaniaGraph(PATH_TO_GLOBAL_GRAPH + 
                                     'interaction_types/combined.ppi.network')

    # Get Global Vertex, Edge, Module Set
    globalVertexSet = getVertexSet(globalGraph)
    globalEdgeSet = getEdgeSet(globalGraph)
    globalModuleSet = getTissueModuleIDSet('global')

    # Get Union Vertex, Edge, Module Set
    unionVertexSet = set()
    unionEdgeSet = set()
    unionModuleSet = set()

    init = True

    # Iterate through tissues
    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        tissueGraph = getTissueSubgraph(tissue)
        vertexSet = getVertexSet(tissueGraph)
        edgeSet = getEdgeSet(tissueGraph)
        moduleSet = getTissueModuleIDSet(tissue)
        
        unionVertexSet = unionVertexSet.union(vertexSet)
        unionEdgeSet = unionEdgeSet.union(edgeSet)
        unionModuleSet = unionModuleSet.union(moduleSet)

    vertexDifference = globalVertexSet.difference(unionVertexSet)
    edgeDifference = globalEdgeSet.difference(unionEdgeSet)
    moduleDifference = globalModuleSet.difference(unionModuleSet)

    f = open(PATH_TO_GLOBAL_LOCAL_DIFF + 'vertex.difference', 'w')
    for vertex in vertexDifference:
        f.write('%s\n' % str(vertex))
        
    f = open(PATH_TO_GLOBAL_LOCAL_DIFF + 'edge.difference', 'w')
    for edge in edgeDifference:
        geneList = list(edge)
        skip = False

        if IGNORE_INHERITED_DIFF:
            for gene in geneList:
                if gene in vertexDifference:
                    skip = True
                    break

            if skip:
                continue

        f.write('%s\t%s\n' % (edge[0], edge[1]))

    f = open(PATH_TO_GLOBAL_LOCAL_DIFF + 'module.difference', 'w')
    for module in moduleDifference:
        geneList = list(module)
        skip = False
 
        if IGNORE_INHERITED_DIFF:
            for gene in geneList:
                if gene in vertexDifference:
                    skip = True
                    break
 
            if skip:
                continue

        f.write(geneList[0])
        for gene in geneList[1:]:
            f.write('\t%s' % gene)
        f.write('\n')
 
        
def correctEdgeDifference():
    f1 = open(PATH_TO_ANALYSIS + 'vertex.difference', 'r')
    vertexSet = set()
    for line in f1:
        gene = line.strip()
        vertexSet.add(gene)
        
    f2 = open(PATH_TO_ANALYSIS + 'edge.difference', 'r')
    for line in f2:
        geneA = edge[0]
        geneB = edge[1]
        if geneA in vertexDifference:
            continue
        if geneB in vertexDifference:
            continue
        f.write('%s\n' % str(edge))

    f = open(PATH_TO_ANALYSIS + 'module.difference', 'r')
    for module in moduleDifference:
        f.write('%s\n' % str(module))
    
def getGlobalLocalVertexAttributes():
    # Open Input File
    inFilename = 'vertex.difference'
    inFilepath = PATH_TO_GLOBAL_LOCAL_DIFF + 'no_double_jeopardy/'
    inFile = open(inFilepath + inFilename, 'r')

    # Open Output File
    outFilename = 'global.local.vertex.attributes'
    outFilepath = PATH_TO_GLOBAL_LOCAL_DIFF
    outFile = open(outFilepath + outFilename, 'w')
    
    globalOnlyVertices = set()

    # Get all vertices that are only global
    for line in inFile:
        gene = parseTabSeparatedLine(line)[0]
        globalOnlyVertices.add(gene)
    
    # Get all vertices for global graph
    globalGraph = getTissueSubgraph('global')
    globalGeneSet = getVertexSet(globalGraph)

    # Write to output file
    for gene in globalGeneSet:
        if gene in globalOnlyVertices:
            outFile.write('%s\t%s\n' % (gene, 'global'))
        else:
            outFile.write('%s\t%s\n' % (gene, 'local'))

    inFile.close()
    outFile.close()
        
    return

def getGlobalLocalEdgeAttributes():
    # Open Input File
    inFilename = 'edge.difference'
    inFilepath = PATH_TO_GLOBAL_LOCAL_DIFF + 'no_double_jeopardy/'
    inFile = open(inFilepath + inFilename, 'r')

    # Open Output File
    outFilename = 'global.local.edge.attributes'
    outFilepath = PATH_TO_GLOBAL_LOCAL_DIFF
    outFile = open(outFilepath + outFilename, 'w')

    globalOnlyEdges = set()

    # Get all edges that are only global
    for line in inFile:
        edgeList = parseTabSeparatedLine(line)
        edgeTuple = tuple(edgeList)
        globalOnlyEdges.add(edgeTuple)

    # Get all edges for global graph
    globalGraph = getTissueSubgraph('global')
    globalEdgeSet = getEdgeSet(globalGraph)

    # Write to output file
    for edge in globalEdgeSet:
        if edge in globalOnlyEdges:
            outFile.write('%s\t%s\t%s\t\n' % (edge[0], 'global', edge[1]))
        else:
            outFile.write('%s\t%s\t%s\t\n' % (edge[0], 'local', edge[1]))
    
    inFile.close()
    outFile.close()

    return

def getGlobalTissueEdgeAttributes(tissue):
    # Open Input File
    inFilename = 'edge.difference'
    inFilepath = PATH_TO_GLOBAL_LOCAL_DIFF + 'no_double_jeopardy/'
    inFile = open(inFilepath + inFilename, 'r')

    # Open Output File
    outFilename = 'global.%s.edge.attributes' % tissue
    outFilepath = PATH_TO_GLOBAL_LOCAL_DIFF
    outFile = open(outFilepath + outFilename, 'w')

    tissueEdges = set()

    # Get all edges that are in tissue
    tissueGraph = getTissueSubgraph(tissue)
    tissueEdgeSet = getEdgeSet(tissueGraph)

    # Get all edges for global graph
    globalGraph = getTissueSubgraph('global')
    globalEdgeSet = getEdgeSet(globalGraph)

    # Write to output file
    for edge in globalEdgeSet:
        if edge in tissueEdgeSet:
            outFile.write('%s\t%s\t%s\t\n' % (edge[0], tissue, edge[1]))
        else:
            outFile.write('%s\t%s\t%s\t\n' % (edge[0], 'global', edge[1]))
    
    inFile.close()
    outFile.close()

    return
