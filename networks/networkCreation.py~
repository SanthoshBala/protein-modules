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


PATH_TO_GENEMANIA_NETWORKS = PATH_TO_NETWORKS + 'genemania/Homo_sapiens.COMBINED/'

PATH_TO_GLOBAL_GRAPH = PATH_TO_GENEMANIA_NETWORKS + 'interaction_types/'
GLOBAL_GRAPH_FILENAME = 'combined.ppi.network'


# - - - - - - - - - - TISSUE SUBGRAPHS - - - - - - - - - - #


# createAllTissueSubgraphs: Creates all tissue subgraphs from global graph,
# including copying the global graph to new location and creating intersection
# graph, composed of proteins that are in all tissue-specific networks.
def createAllTissueSubgraphs():
    # Create Global Graph
    createGlobalSubgraph()
    if PRINT_PROGRESS:
        print 'global'
    
    # Create Intersection Graph
    createIntersectionSubgraph()
    if PRINT_PROGRESS:
        print 'intersection'

    # Print Tissue Subgraphs
    for tissue in FUNCTIONAL_TISSUE_LIST:
        createTissueSubgraph(tissue)    
        if PRINT_PROGRESS:
            print tissue
    
    return

# createGlobalSubgraph: Automates creating copy of global graph in tissue
# subgraph directory to simplify downstream iteration through tissues.
def createGlobalSubgraph():
    oldFilePath = PATH_TO_GLOBAL_GRAPH + GLOBAL_GRAPH_FILENAME
    newFileName = 'global.' + GLOBAL_GRAPH_FILENAME
    newFilePath = PATH_TO_TISSUE_SUBGRAPHS + newFileName

    shutil.copy(oldFilePath, newFilePath)
    
    return

# createIntersectionSubgraph: Creates graph that represents intersection
# of all tissue-specific subgraphs.
def createIntersectionSubgraph():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Input File
    inFile = open(PATH_TO_GLOBAL_GRAPH + GLOBAL_GRAPH_FILENAME, 'r')

    numTissues = len(FUNCTIONAL_TISSUE_LIST)

    headerLine = True
    
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
        recordA = ngtmDB.find_one( { 'gene_id' : geneA } )
        recordB = ngtmDB.find_one( { 'gene_id' : geneB } )
        if not recordA or not recordB:
            continue

        tissuesA = recordA.get('tissue_list')
        tissuesB = recordB.get('tissue_list')
        if len(tissuesA) < numTissues:
            continue
        if len(tissuesB) < numTissues:
            continue

        outFile.write(line)     

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return

# createTissueSubgraph: Creates subgraph for <tissue> from global graph.
def createTissueSubgraph(tissue):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ntgmDB = db.normTissueGeneMap

    # Open Input File
    inFile = open(PATH_TO_GLOBAL_GRAPH + GLOBAL_GRAPH_FILENAME, 'r')
    outFile = open(PATH_TO_TISSUE_SUBGRAPHS + tissue + '.' + 
                   GLOBAL_GRAPH_FILENAME)

    tissueRecord = ntgmDB.find_one( { 'tissue' : tissue } )
    geneDict = tissueRecord.get(DESIRED_NOMENCLATURE)
    
    headerLine = True

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
        
        # If 1 of genes not in tissue, continue
        if geneDict.get(geneA) < ARRAY_EXPRESSION_THRESHOLD:
            continue
        if geneDict.get(geneB) < ARRAY_EXPRESSION_THRESHOLD:
            continue

        outFile.write(line)

    # Close Files and DB Connection
    inFile.close()
    outFile.close()
    cli.close()

    return
