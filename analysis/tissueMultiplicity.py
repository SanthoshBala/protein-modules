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

# Utility Imports
from common.matrices import *


# - - - - - - - - - - - TISSUE SPECIFICITY - - - - - - - - - - #

PATH_TO_TISSUE_MULTIPLICITY = PATH_TO_ANALYSIS + 'tissue_multiplicity/' 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# getProteinTissueMultiplicity: Write file containing names of proteins
# and number of tissues in which each one is found.
def getProteinTissueMultiplicity():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Output File
    outFilePath = PATH_TO_TISSUE_MULTIPLICITY
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

# getInteractionTissueMultiplicity: Write file containing two dimensional
# histogram, showing prevalence of protein interaction pairs.
def getInteractionTissueMultiplicity():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Output File
    outFilePath = PATH_TO_TISSUE_MULTIPLICITY
    outFileName = 'interaction.tissue.multiplicity'
    outFile = open(outFilePath + outFileName, 'w')

    # Get Matrix
    numTissues = len(FUNCTIONAL_TISSUE_LIST)
    matrix = constructSquareIntegerMatrix(numTissues + 1)
    
    # Get Global Interaction Set
    globalGraph = getTissueSubgraph('global')
    interactionSet = getInteractionSet(globalGraph)

    # Iterate through Interactions
    for interaction in interactionSet:
        if PRINT_PROGRESS:
            print interaction

        proteinA = interaction[0]
        proteinB = interaction[0]

        recordA = ngtmDB.find_one( { 'gene_id' : proteinA } )
        if not recordA:
            continue
        recordB = ngtmDB.find_one( { 'gene_id' : proteinB } )
        if not recordB:
            continue

        tissuesA = len(recordA.get('tissue_list'))
        tissuesB = len(recordB.get('tissue_list'))

        matrix[tissuesA][tissuesB] = matrix[tissuesA][tissuesB] + 1
        matrix[tissuesB][tissuesA] = matrix[tissuesB][tissuesA] + 1

    # Write Matrix to File
    for i in range(numTissues + 1):
        for j in range(numTissues + 1):
            outFile.write('%d\t%d\t%d\n' % (i, j, matrix[i][j]))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return

# getModuleTissueMultiplicity
def getModuleTissueMultiplicity():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilePath = PATH_TO_TISSUE_MULTIPLICITY
    outFileName = 'module.tissue.multiplicity'
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through <modDB>
    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        tissueList = moduleRecord.get('tissue_list')
        tissueMultiplicity = len(tissueList)
        
        if 'global' in tissueList:
            if len(tissueList) == 1:
                continue
            else:
                tissueMultiplicity = tissueMultiplicity - 1
        if 'intersection' in tissueList:
            tissueMultiplicity = len(FUNCTIONAL_TISSUE_LIST)
            
        outFile.write('%s\t%d\n' % (moduleID, tissueMultiplicity))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return

# getModuleProteinUniversality
def getModuleProteinUniversality():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilePath = PATH_TO_TISSUE_MULTIPLICITY
    outFileName = 'module.protein.universality'
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through <modDB>
    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        universality = moduleRecord.get('protein_universality')
            
        outFile.write('%s\t%f\n' % (moduleID, universality))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return

# getModuleHousekeepingIndex
def getModuleHousekeepingIndex():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilePath = PATH_TO_TISSUE_MULTIPLICITY
    outFileName = 'module.housekeeping.index'
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through <modDB>
    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        housekeeping = moduleRecord.get('housekeeping_index')
            
        outFile.write('%s\t%f\n' % (moduleID, housekeeping))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
