#! /usr/bin/python

import settings
from settings import *
import pymongo
from pymongo import MongoClient
import os
from graphIO import *
from moduleUtil import *
from moduleTopology import collapseModule

PATH_TO_FILE = '/home/santhosh/workspace/thesis/data/networks/genemania/Homo_sapiens.COMBINED/'
FILENAME = 'COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt'

OUT_BASE_FILENAME = '%s.tissueSubgraph.analysis'

# - - - - - GRAPH GENERATION - - - - - #

# generateAllCollapsedTissueSubgraphs():
def generateAllCollapsedTissueSubgraphs():
    # Use Augmented Tissue List
    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    # Iterate through all Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        generateCollapsedTissueSubgraph(tissue)
    
# generateCollapsedTissueSubgraph
def generateCollapsedTissueSubgraph(tissue):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    moduleDB = db.modules

    # Get Tissue Subgraph
    graph = getTissueSubgraph(tissue)

    # Get Module Set
    moduleSet = getTissueModuleGeneSet(tissue)

    # Iterate through Modules for Tissue and Collapse within Graph
    for module in moduleSet:
        geneList = list(module)
        collapseModule(graph, geneList)

    # Write Collapsed Graph to Disk
    outFilename = GENEMANIA_COLLAPSED_TISSUE_BASE_FILENAME % \
        (tissue, DESIRED_INTERACTION_TYPE)
    outFilePath = PATH_TO_COLLAPSED_TISSUE_SUBGRAPHS + outFilename
    writeGenemaniaGraphToDisk(graph, outFilePath)

    cli.close()
    return

# - - - - - GRAPH ANALYSIS - - - - - #

# analyzeAllTissueSubgraphs
def analyzeAllTissueSubgraphs():

    for interactionType in GENEMANIA_INTERACTION_TYPES:
        print interactionType
        
        outputFilename = OUT_BASE_FILENAME % interactionType
        outputFile = open(PATH_TO_ANALYSIS + outputFilename, 'w')
        outputFile.write('Analyzing all tissue subgraphs...\n\n')
        outputFile.write('Tissue\tV\tE\tDeg\tDensity\n')
        
    # Analyze all tissues for this interaction type
        for tissue in CANONICAL_TISSUE_LIST:
            inputFilename = GENEMANIA_GLOBAL_BASE_FILENAME % interactionType
            generateTissueSubgraph(tissue, 'primary_gene_id', inputFilename)
            numV, numE, avgD, density = analyzeTissueGraph(tissue + '.' + inputFilename)
            outputFile.write('%s\t%d\t%d\t%f\t%f\n' % \
                                 (tissue, numV, numE, avgD, density))

            outputFile.flush()
            print '\t%s' % tissue

# analyzeTissueGraph
def analyzeTissueGraph(filename):
    f = open(PATH_TO_TISSUE_SUBGRAPHS + filename, 'r')

    firstLine = True

    adjList = {}
    
    for line in f:
        if firstLine:
            firstLine = False
            continue

        lineTabFields = line.split('\t')
        lineFields = [s.strip() for s in lineTabFields]

        geneA = lineFields[0]
        geneB = lineFields[1]
        
        if adjList.get(geneA):
            adjList[geneA].append(geneB)
        else:
            adjList.update( { geneA : [geneB] } )
            
        if adjList.get(geneB):
            adjList[geneB].append(geneA)
        else:
            adjList.update( { geneB : [geneA] } )

    # Iterate through adjacency list and print output
    numNodes = len(adjList.keys())
    numEdges = 0
    for key, val in adjList.iteritems():
        numEdges = numEdges + len(val)
        
    numEdges = numEdges / 2

    avgDegree = float(2*numEdges)/float(numNodes)
    density = float(2*numEdges)/float((numNodes*(numNodes-1)))

    f.close()
    
    return numNodes, numEdges, avgDegree, density
    
        
