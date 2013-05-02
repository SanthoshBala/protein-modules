#! /usr/bin/python

# proteomeCoverage.py
# Author: Santhosh Balasubramanian
# Created: April 14, 2013
# Last Modified: April 14, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphIO import *
from graphs.graphUtil import *

# Utility Imports
from common.strings import *
from common.statistics import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

PATH_TO_PROTEOME_COVERAGE = PATH_TO_ANALYSIS + 'proteome_coverage/'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# getProteomeCoverage: Creates file containing number of proteins per tissue
# according to GeneMANIA and Su et al. data.
def getProteomeCoverage(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        ngtmDB = db.shuffleNormGeneTissueMap
    else:
        ngtmDB = db.normGeneTissueMap

    # Open Output File
    if shuffle:
        outFileName = 'shuffle.proteome.coverage'
        outFilePath = PATH_TO_PROTEOME_COVERAGE
    else:
        outFileName = 'proteome.coverage'
        outFilePath = PATH_TO_PROTEOME_COVERAGE
    outFile = open(outFilePath + outFileName, 'w')

    outFile.write('Tissue\tMicroarray\tPPI\tCoverage\n')

    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue

        # Get Num Proteins According to Microarray
        if tissue == 'global':
            numArrayProteins = ngtmDB.count()
        elif tissue == 'intersection':
            numTissues = len(FUNCTIONAL_TISSUE_LIST)
            numArrayProteins = ngtmDB.find( { 'tissue_list' : 
                                              { '$size' : 
                                                numTissues } } ).count()
        else:
            numArrayProteins = ngtmDB.find( { 'tissue_list' : tissue } ).count()

        # Get Num Proteins According to PPI
        tissueSubgraph = getTissueSubgraph(tissue)
        numPPIProteins = tissueSubgraph.num_vertices()
        
        coverage = float(numPPIProteins)/float(numArrayProteins)

        outFile.write('%s\t%d\t%d\t%f\n' % (tissue, numArrayProteins, 
                                        numPPIProteins, coverage))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return


# computeCoverageStatistics: Computes mean coverage and standard deviation.
def computeCoverageStatistics(shuffle = False):
    
    if shuffle:
        inFileName = 'shuffle.proteome.coverage'
        outFileName = 'shuffle.proteome.coverage.statistics'
    else:
        inFileName = 'proteome.coverage'
        outFileName = 'proteome.coverage.statistics'

    # Open Input File
    inFilePath = PATH_TO_PROTEOME_COVERAGE
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    outFilePath = PATH_TO_PROTEOME_COVERAGE
    outFile = open(outFilePath + outFileName, 'w')


    coverageList = list()
    headerLine = True

    # Iterate through <inFile>
    for line in inFile:
        if headerLine:
            headerLine = False
            continue

        lineFields = parseTabSeparatedLine(line)
        coverage = lineFields[3]
        coverageList.append(float(coverage))
    
    # Get Mean and Std Dev
    avg = mean(coverageList)[0]
    dev = stdDev(coverageList)

    outFile.write('Mean Coverage\t%f\n' % avg)
    outFile.write('Standard Deviation\t%f\n' % dev)

    # Close Files
    inFile.close()
    outFile.close()

    return

    
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
    outFilePath = PATH_TO_PROTEOME_COVERAGE
    outFile = open(outFilePath + outFileName, 'w')

    outFile.write('Property\tGlobal\tLocal\tDifference\n')
    outFile.write('|V|\t%d\t%d\t%f\n' % (globalV, unionV, diffV))
    outFile.write('|E|\t%d\t%d\t%f\n' % (globalE, unionE, diffE))
    outFile.write('|M|\t%d\t%d\t%f\n' % (globalM, unionM, diffM))
    
    # Close File and DB Connection
    outFile.close()    
    cli.close()

    return
