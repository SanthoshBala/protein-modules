#! /usr/bin/python

# tissueSpecificity.py
# Author: Santhosh Balasubramanian
# Created: February 27, 2013
# Last Modified: April 14, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphIO import *
from graphs.graphUtil import *

# - - - - - - - - - - - TISSUE SPECIFICITY - - - - - - - - - - #

PATH_TO_TISSUE_SPECIFICITY = PATH_TO_ANALYSIS + 'tissue_specificity/' 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# compareGlobalAndLocalNetworks: Computes difference between global and union
# of tissue networks in terms of vertex, edge, and module composition.
def compareGlobalAndLocalNetworks():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Get Global Graph
    globalGraph = getTissueSubgraph('global')

    # Get Tissue Union Graph
    unionGraph = getTissueUnionGraph()

    # Get Vertex Sets
    globalVertexSet = getGeneSet(globalGraph)
    unionVertexSet = getGeneSet(unionGraph)
    diffVertexSet = globalVertexSet.difference(unionVertexSet)
    globalV = len(globalVertexSet)
    unionV = len(unionVertexSet)

    # Get Edge Sets, Avoiding Double Jeopardy
    globalEdgeSet = getInteractionSet(globalGraph)
    unionEdgeSet = getInteractionSet(unionGraph)
    diffEdgeSet = globalEdgeSet.difference(unionEdgeSet)
    badDiffEdgeSet = set()
    for edge in diffEdgeSet:
        if edge[0] in diffVertexSet:
            badDiffEdgeSet.add(edge)      
        if edge[1] in diffVertexSet:
            badDiffEdgeSet.add(edge)
       
    globalE = len(globalEdgeSet) - len(badDiffEdgeSet)
    unionE = len(unionEdgeSet)

    # Get Module Sets, Avoiding Double Jeopardy
    globalModuleList = modDB.find( { 'tissue_list' : [ 'global' ] } )
    badModuleSet = set()
    for moduleRecord in globalModuleList:
        # Iterate through Genes
        for gene in moduleRecord.get('gene_list'):
            if gene in diffVertexSet:
                badModuleSet.add(moduleRecord.get('module_id'))
                break

    globalM = modDB.find().count() - len(badModuleSet)
    unionM = modDB.find().count() - globalModuleList.count()

    diffV = float(unionV)/float(globalV)
    diffE = float(unionE)/float(globalE)
    diffM = float(unionM)/float(globalM)

    # Write Output
    outFileName = 'global.local.difference'
    outFilePath = PATH_TO_TISSUE_SPECIFICITY
    outFile = open(outFilePath + outFileName, 'w')

    outFile.write('Property\tGlobal\tLocal\tDifference\n')
    outFile.write('|V|\t%d\t%d\t%f\n' % (globalV, unionV, diffV))
    outFile.write('|E|\t%d\t%d\t%f\n' % (globalE, unionE, diffE))
    outFile.write('|M|\t%d\t%d\t%f\n' % (globalM, unionM, diffM))
    
    # Close File and DB Connection
    outFile.close()    
    cli.close()

    return


# getProteinTissueMultiplicity: Write file containing names of proteins
# and number of tissues in which each one is found.
def getProteinTissueMultiplicity():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Output File
    outFilePath = PATH_TO_TISSUE_SPECIFICITY
    outFileName = 'protein.tissue.multiplicity'
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through DB
    for proteinRecord in ngtmDB.find():
        proteinID = proteinRecord.get('gene_id')
        tissueMultiplicity = len(proteinRecord.get('tissue_list'))
        outFile.write('%s\t%d\n' % (proteinID, tissueMultiplicity))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
