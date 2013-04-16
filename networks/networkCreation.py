#! /usr/bin/python

# networkCreation.py
# Author: Santhosh Balasubramanian
# Created: April 10, 2013
# Last Modified: April 10, 2013


# Python Imports
import os
import shutil

# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.strings import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

PATH_TO_GLOBAL_GRAPH = PATH_TO_GENEMANIA_NETWORKS + 'interaction_types/'
GLOBAL_GRAPH_FILENAME = 'combined.ppi.network'

# - - - - - - - - - - TISSUE SUBGRAPHS - - - - - - - - - - #


# createGlobalNetwork: Automates creating copy of global graph in tissue
# subgraph directory to simplify downstream iteration through tissues.
# If <shuffle>, uses the randomized tissue expression.
def createGlobalNetwork(shuffle = False):
    oldFilePath = PATH_TO_GLOBAL_GRAPH + GLOBAL_GRAPH_FILENAME

    if shuffle:
        newFileName = SHUFFLE_TISSUE_SUBGRAPH_BASE_FILENAME % 'global'
        newFilePath = PATH_TO_SHUFFLE_TISSUE_SUBGRAPHS + newFileName
    else:
        newFileName = 'global.' + GLOBAL_GRAPH_FILENAME
        newFilePath = PATH_TO_TISSUE_SUBGRAPHS + newFileName

    shutil.copy(oldFilePath, newFilePath)

    if PRINT_PROGRESS:
        print 'global'

    return

# createIntersectionNetwork: Creates graph that represents intersection
# of all tissue-specific subgraphs. If <shuffle>, uses randomized tissue
# expression.
def createIntersectionNetwork(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        gtmDB = db.shuffleNormGeneTissueMap
    else:
        gtmDB = db.normGeneTissueMap

    # Open Input File
    inFile = open(PATH_TO_GLOBAL_GRAPH + GLOBAL_GRAPH_FILENAME, 'r')

    # Open Output File
    if shuffle:
        outFilePath = PATH_TO_SHUFFLE_TISSUE_SUBGRAPHS
        outFileName = SHUFFLE_TISSUE_SUBGRAPH_BASE_FILENAME % ('intersection')
    else:
        outFilePath = PATH_TO_TISSUE_SUBGRAPHS
        outFileName = GENEMANIA_TISSUE_BASE_FILENAME % ('intersection', 'ppi')
    outFile = open(outFilePath + outFileName, 'w')
        
    numTissues = len(FUNCTIONAL_TISSUE_LIST)

    headerLine = True
    
    # Get Housekeeping Proteins
    numTissues = len(FUNCTIONAL_TISSUE_LIST)
    housekeepingSet = set()
    housekeepingRecords = gtmDB.find( { 'tissue_list' : 
                                        { '$size' : numTissues } }, 
                                      { 'module_id' : 1 } )
    for record in housekeepingRecords:
        housekeepingSet.add(record.get('module_id'))

    # Iterate through <inFile>
    for line in inFile:
        if headerLine:
            headerLine = False
            outFile.write(line)
            continue

        lineTabFields = parseTabSeparatedLine(line)

        geneA = lineTabFields[0]
        geneB = lineTabFields[1]
        confidence = lineTabFields[2]
        
        # If 1 of genes not in all tissues, continue
        if geneA not in housekeepingSet:
            continue
        if geneB not in housekeepingSet:
            continue

        outFile.write(line)     

    # Close File and DB Connection
    inFile.close()
    outFile.close()
    cli.close()

    return

# createTissueNetworks: Creates all tissue subgraphs from global graph,
# including copying the global graph to new location and creating intersection
# graph, composed of proteins that are in all tissue-specific networks.
def createTissueNetworks(minimum = 0, maximum = 80, shuffle = False):
    # Create Tissue Subgraphs
    for tissue in FUNCTIONAL_TISSUE_LIST[minimum:maximum]:
        createTissueNetwork(tissue, shuffle = shuffle)    
        if PRINT_PROGRESS:
            print tissue
    
    return

# createTissueNetwork
def createTissueNetwork(tissue, shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        ntgmDB = db.shuffleNormTissueGeneMap
    else:
        ntgmDB = db.normTissueGeneMap


    # Open Input File
    inFile = open(PATH_TO_GLOBAL_GRAPH + GLOBAL_GRAPH_FILENAME, 'r')
    
    # Open Output File
    if shuffle:
        outFilePath = PATH_TO_SHUFFLE_TISSUE_SUBGRAPHS
        outFileName = SHUFFLE_TISSUE_SUBGRAPH_BASE_FILENAME % (tissue)
    else:
        outFilePath = PATH_TO_TISSUE_SUBGRAPHS
        outFileName = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, 'ppi')
    outFile = open(outFilePath + outFileName, 'w')

    headerLine = True

    tissueRecord = ntgmDB.find_one( { 'tissue' : tissue } )
    geneDict = tissueRecord.get('primary_gene_id')

    # Iterate through <inFile>
    for line in inFile:
        if headerLine:
            headerLine = False
            outFile.write(line)
            continue

        lineTabFields = parseTabSeparatedLine(line)

        geneA = lineTabFields[0]
        geneB = lineTabFields[1]
        
        # Get DB Records
        expressionA = geneDict.get(geneA)
        if expressionA <= float(ARRAY_EXPRESSION_THRESHOLD):
            continue
        expressionB = geneDict.get(geneB)
        if expressionB <= float(ARRAY_EXPRESSION_THRESHOLD):
            continue

        outFile.write(line)

    # Close Files and DB Connection
    inFile.close()
    outFile.close()
    cli.close()

    return
    


# - - - - - - - - - - COLLAPSED TISSUE SUBGRAPHS - - - - - - - - - - #


# createAllCollapsedTissueSubgraphs: Creates collapsed graphs for all tissues,
# including global graph and intersection graph.
def createAllCollapsedTissueSubgraphs():
    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue
        createCollapsedTissueSubgraph(tissue)

    return


# createCollapsedTissueSubgraph: Collapses all genes in a module to a single
# node, and writes resulting graph to disk.
def createCollapsedTissueSubgraph(tissue):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    moduleDB = db.modules

    # Get Tissue Subgraph
    graph = getTissueSubgraph(tissue)

    # Get Module Set
    moduleSet = getTissueRawModuleSet

    # Iterate through Modules for Tissue and Collapse
    for module in moduleSet:
        geneList = list(module)
        collapseModule(graph, geneList)

    # Write Collapsed Graph to Disk
    outFileName = GENEMANIA_COLLAPSED_TISSUE_BASE_FILENAME % (tissue, 'ppi')
    outFilePath = PATH_TO_COLLAPSED_TISSUE_SUBGRAPHS + outFileName
    writeGenemaniaGraphToDisk(graph, outFilePath)

    # Close DB Connection
    cli.close()

    return



