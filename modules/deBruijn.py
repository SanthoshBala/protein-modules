#! /usr/bin/python

# deBruijn.py
# Author: Santhosh Balasubramanian
# Created: March 16, 2013
# Last Modified: March 27, 2013


# Python Imports
from copy import *

# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.matrices import *


# - - - - - - - - - - ALGORITHM - - - - - - - - - - #


# De Bruijn graphs are a particular graph construction used for modeling
# the relationships between sequences. The code here is used for constructing
# modified De Bruijn graphs that can compare sets. In these graphs, sets
# are represented by nodes, and two nodes (sets) are adjacent if they
# are similar in some way. The three similarity types are exchange, addition,
# and removal. 
#
# A is exchange similar to B iff |A| = |B| and A and B differ by 1 element.
# A is addition similar to B iff |A| = |B| + 1 and B is a subset of A.
# A is removal similar to B iff |A| = |B| - 1 and A is a subset of B.
#
# Note that these relationships are asymmetric, and that the graph is
# therefore a directed graph. The functions here allow for the simple
# construction of De Bruijn set graphs for comparing PPI modules.


# - - - - - - - - - - SIMILAR MODULE RETRIEVAL - - - - - - - - - - #


# getDeBruijnSimilarModules: Retrieves all modules that are "similar"
# to <moduleID> in a dynamic De Bruijn sense. That is, they are "similar"
# if the two modules can be related by a single exchange, addition, or
# removal operation.
def getDeBruijnSimilarModules(moduleID):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    deBruijnAdjacentSet = set()

    # Get Similar Modules by Type
    exchangeMods = getExchangeSimilarModules(moduleID)
    additionMods = getAdditionSimilarModules(moduleID)
    removalMods = getRemovalSimilarModules(moduleID)

    # Get Union of Similar Modules
    deBruijnAdjacentSet = deBruijnAdjacentSet.union(exchangeMods)
    deBruijnAdjacentSet = deBruijnAdjacentSet.union(additionMods)
    deBruijnAdjacentSet = deBruijnAdjacentSet.union(removalMods)

    # Close DB Connection
    cli.close()

    return deBruijnAdjacentSet

    
# getExchangeSimilarModules: Get all modules that are found in any
# tissue of <tissueSet> and are "exchange similar" to <moduleID>.
# Two modules are exchange similar if they have the same size
# and differ by only one element.
def getExchangeSimilarModules(moduleID, tissueList = AUGMENTED_TISSUE_LIST):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    similarMods = set()
    tissueSet = set(tissueList)

    # Get Gene List for <moduleID>
    moduleRecord = modDB.find_one( { 'module_id' : moduleID } )
    moduleGeneList = moduleRecord.get('gene_list')
    numGenes = len(moduleGeneList)

    # Iterate through Genes in <moduleID>
    for i in range(numGenes):
        truncGeneList = moduleGeneList[0:i] + moduleGeneList[i+1:numGenes]

        # Substitution Similar Modules Have Same # of Genes as <moduleID>
        substitutionMods = modDB.find( {'gene_list':{'$all':truncGeneList}} )

        # Iterate through Possible Similar Modules
        for module in substitutionMods:
            # Ignore <moduleID>
            if module.get('module_id') == moduleID:
                continue

            # Check <module> in Desired Tissues
            moduleTissueSet = set(module.get('tissue_list'))
            if len(moduleTissueSet.intersection(tissueSet)) == 0:
                continue

            # Check Module Size
            if len(module.get('gene_list')) == numGenes:
                similarMods.add(module.get('module_id'))

    # Close DB Connection
    cli.close()

    return similarMods

# getAdditionSimilarModules: Get all modules that are found in any
# tissue of <tissueSet> and are "addition similar" to <moduleID>.
# A is addition similar to B if |A| = |B| + 1, and B is a subset of A.
def getAdditionSimilarModules(moduleID, tissueList = AUGMENTED_TISSUE_LIST):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    additionMods = set()
    tissueSet = set(tissueList)

    # Get Gene List for "moduleID"
    moduleRecord = modDB.find_one( { 'module_id' : moduleID } )
    moduleGeneList = moduleRecord.get('gene_list')
    numGenes = len(moduleGeneList)

    # Insertion Similar Have 1 More Gene than <moduleID>
    similarMods = modDB.find( { 'gene_list' : { '$all': moduleGeneList } } )

    # Iterate through Possible Similar Modules
    for module in similarMods:
        # Ignore "moduleID" itself
        if module.get('module_id') == moduleID:
            continue
        
        # Check <module> in Desired Tissues
        moduleTissueSet = set(module.get('tissue_list'))
        if len(moduleTissueSet.intersection(tissueSet)) == 0:
            continue

        # Check Module Size
        if len(module.get('gene_list')) == numGenes + 1:
            additionMods.add(module.get('module_id'))
            
    # Close DB Connection
    cli.close()

    return additionMods

# getRemovalSimilarModules: Get all modules that are found in any
# tissue of <tissueSet> and are "removal similar" to <moduleID>.
# A is removal similar to B if |A| = |B| - 1, and A is a subset of B.
def getRemovalSimilarModules(moduleID, tissueList = AUGMENTED_TISSUE_LIST):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    removalMods = set()
    tissueSet = set(tissueList)

    # Get Gene List for <moduleID>
    moduleRecord = modDB.find_one( { 'module_id' : moduleID } )
    moduleGeneList = moduleRecord.get('gene_list')
    numGenes = len(moduleGeneList)
    
    # Iterate through Genes in <moduleID>
    for i in range(numGenes):
        truncGeneList = moduleGeneList[0:i] + moduleGeneList[i+1:numGenes]

        # Deletion Similar modules Have 1 less Gene as <moduleID>
        similarMods = modDB.find( {'gene_list':{ '$all' : truncGeneList } } )

        # Iterate through Possible Similar Modules
        for module in similarMods:
            # Ignore <moduleID> itself
            if module.get('module_id') == moduleID:
                continue
            
            # Check <module> in Desired Tissues
            moduleTissueSet = set(module.get('tissue_list'))
            if len(moduleTissueSet.intersection(tissueSet)) == 0:
                continue
            
            # Check Module Size
            if len(module.get('gene_list')) == numGenes - 1:
                removalMods.add(module.get('module_id'))

    # Close DB Connection
    cli.close()

    return removalMods


# - - - - - - - - - - DE BRUIJN SET GRAPH CREATION - - - - - - - - - - #


# createDeBruijnSetGraph: Creates graph in which nodes are modules, and edges
# represent similarities between modules. Each edge is annotated with the
# similarity type. The graph is written to <outFileName>.
def createDeBruijnSetGraph(outFileName, tissueList = AUGMENTED_TISSUE_LIST):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    tissueSet = set(tissueList)

    # Open Output File
    outFilePath = PATH_TO_MODULE_SIMILARITY + outFileName
    outFile = open(outFilePath, 'w')

    # Iterate through Modules
    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        moduleTissueSet = set(moduleRecord.get('tissue_list'))

        if PRINT_PROGRESS:
            print moduleID

        # Check Module in Tissue Set
        if len(tissueSet.intersection(moduleTissueSet)) == 0:
            continue

        # Check Module Size
        if len(moduleRecord.get('gene_list')) < MODULE_SIZE_THRESHOLD:
            continue

        # Get Similar Modules for <moduleID>
        exchangeSimMods = getExchangeSimilarModules(moduleID, tissueList)
        additionSimMods = getAdditionSimilarModules(moduleID, tissueList)
        removalSimMods = set()

        # Ignore Deletion Similar Modules for 4-Protein Modules
        if len(moduleRecord.get('gene_list')) != MODULE_SIZE_THRESHOLD:
            removalSimMods = getRemovalSimilarModules(moduleID, tissueSet)

        # Write Edges to File
        for exchange in exchangeSimMods:
            outFile.write('%s\t%s\texchange\n' % (moduleID, exchange))
        for addition in additionSimMods:
            outFile.write('%s\t%s\taddition\n' % (moduleID, addition))
        for removal in removalSimMods:
            outFile.write('%s\t%s\tremoval\n' % (moduleID, removal))

        # Write Singleton Nodes with Self Edges
        if len(substitutionSimMods) == 0:
            if len(insertionSimMods) == 0:
                if len(deletionSimMods) == 0:
                    outFile.write('%s\t%s\tself\n' % (moduleID, moduleID))

    # Close File and DB Connection
    cli.close()
    outFile.close()

    return


# - - - - - - - - - - DE BRUIJN SET GRAPH ANALYSIS - - - - - - - - - - #


# createEmbryonicDeBruijnAdditionMatrix: Creates matrix representing
# the number of addition similarities between vertices of two types. The
# "type" is the germ layer hash of the modules representing the two nodes.
# i.e. (MESO)-->(MESO+ENDO) increments M[2][6] by 1.
def createEmbryonicDeBruijnAdditionMatrix(inFileName):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Input File
    inFilePath = PATH_TO_MODULE_SIMILARITY + inFileName
    inFile = open(inFilePath, 'r')

    additionMatrix = constructSquareIntegerMatrix(8)

    # Iterate through <inFile>
    for line in inFile:
        # Ignore non-Addition Edges
        lineFields = parseTabSeparatedLine(line)
        if lineFields[2] != 'addition':
            continue

        # Get Germ Layer Hashes
        hashA = getGermLayerHash(lineFields[0])
        hashB = getGermLayerHash(lineFields[1])

        additionMatrix[hashA][hashB] = additionMatrix[hashA][hashB] + 1

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return

# createEmbryonicDeBruijnRemovalMatrix: Creates matrix representing
# the number of removal similarities between vertices of two types. The
# "type" is the germ layer hash of the modules representing the two nodes.
# i.e. (MESO)-->(MESO+ENDO) increments M[2][6] by 1.
def constructEmbryonicDeBruijnRemovalMatrix(inFileName):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Input File
    inFilePath = PATH_TO_MODULE_SIMILARITY + inFileName
    inFile = open(inFilePath, 'r')

    removalMatrix = constructSquareIntegerMatrix(8)

    # Iterate through <inFile>
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        # Ignore non-Removal Edges
        if lineFields[2] != 'removal':
            continue

        # Get Germ Layer Hashes
        hashA = getGermLayerHash(lineFields[0])
        hashB = getGermLayerHash(lineFields[1])

        removalMatrix[hashA][hashB] = removalMatrix[hashA][hashB] + 1


    # Close File and DB Connection
    inFile.close()
    cli.close()

    return
    
# createEmbryonicDeBruijnExchangeMatrix: Creates matrix representing
# the number of exchange similarities between vertices of two types. The
# "type" is the germ layer hash of the modules representing the two nodes.
# i.e. (MESO)-->(MESO+ENDO) increments M[2][6] by 1.
def constructEmbryonicDeBruijnExchangeMatrix(inFileName):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Input File
    inFilePath = PATH_TO_MODULE_SIMILARITY + inFileName
    inFile = open(inFilePath, 'r')

    exchangeMatrix = constructSquareIntegerMatrix(8)

    # Iterate through <inFile>
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        # Ignore non-Exchange Edges
        if lineFields[2] != 'exchange':
            continue

        # Get Germ Layer Hashes
        hashA = getGermLayerHash(lineFields[0])
        hashB = getGermLayerHash(lineFields[1])

        exchangeMatrix[hashA][hashB] = exchangeMatrix[hashA][hashB] + 1

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return

# getDeBruijnEdgeTypeCounts: Returns # of Exchanges, Additions, and Removals.
def getDeBruijnEdgeTypeCounts(inFileName):
    # Open Input File
    inFilePath = PATH_TO_MODULE_SIMILARITY + inFileName
    inFile = open(inFilePath, 'r')

    numExchanges = 0
    numAdditions = 0
    numRemovals = 0

    # Iterate through Graph
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        if lineFields[2] == 'addition':
            numInsertions = numAdditions + 1
        if lineFields[2] == 'removal':
            numDeletions = numRemovals + 1
        if lineFields[2] == 'exchange':
            numSubstitutions = numExchanges + 1

    return numExchanges, numAdditions, numRemovals
