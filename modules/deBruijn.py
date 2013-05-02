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

# Module Imports
from modules.moduleUtil import *

# Utility Imports
from common.strings import *
from common.matrices import *

# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

PATH_TO_DEBRUIJN_GRAPH = PATH_TO_ANALYSIS + 'debruijn/'
DEBRUIJN_GRAPH_FILENAME = 'debruijn.similarity.graph'

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
def getExchangeSimilarModules(moduleID, tissueList = AUGMENTED_TISSUE_LIST, shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
    else:
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

        #  Similar Modules Have Same # of Genes as <moduleID>
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
def getAdditionSimilarModules(moduleID, tissueList = AUGMENTED_TISSUE_LIST, shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
    else:
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
def getRemovalSimilarModules(moduleID, tissueList = AUGMENTED_TISSUE_LIST, shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
    else:
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
def createDeBruijnSetGraph(tissueList = AUGMENTED_TISSUE_LIST, shuffleVal = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffleVal:
        modDB = db.shuffleModules
    else:
        modDB = db.modules

    tissueSet = set(tissueList)

    # Open Output File
    if shuffleVal:
        outFileName = 'shuffle.' + DEBRUIJN_GRAPH_FILENAME
    else:
        outFileName = DEBRUIJN_GRAPH_FILENAME
    outFilePath = PATH_TO_DEBRUIJN_GRAPH
    outFile = open(outFilePath + outFileName, 'w')

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
        exchangeSimMods = getExchangeSimilarModules(moduleID, tissueList, shuffle = shuffleVal)
        additionSimMods = getAdditionSimilarModules(moduleID, tissueList, shuffle = shuffleVal)
        removalSimMods = set()

        # Ignore Deletion Similar Modules for 4-Protein Modules
        if len(moduleRecord.get('gene_list')) != MODULE_SIZE_THRESHOLD:
            removalSimMods = getRemovalSimilarModules(moduleID, tissueSet, shuffle = shuffleVal)

        # Write Edges to File
        for exchange in exchangeSimMods:
            outFile.write('%s\t%s\texchange\n' % (moduleID, exchange))
        for addition in additionSimMods:
            outFile.write('%s\t%s\taddition\n' % (moduleID, addition))
        for removal in removalSimMods:
            outFile.write('%s\t%s\tremoval\n' % (moduleID, removal))

        # Write Singleton Nodes with Self Edges
        if len(exchangeSimMods) == 0:
            if len(additionSimMods) == 0:
                if len(removalSimMods) == 0:
                    outFile.write('%s\t%s\tself\n' % (moduleID, moduleID))

    # Close File and DB Connection
    cli.close()
    outFile.close()

    return


# - - - - - - - - - - DE BRUIJN SET GRAPH ANALYSIS - - - - - - - - - - #


# getEmbryonicDeBruijnAdditionMatrix: Creates matrix representing
# the number of addition similarities between vertices of two types. The
# "type" is the germ layer hash of the modules representing the two nodes.
# i.e. (MESO)-->(MESO+ENDO) increments M[2][6] by 1.
def getEmbryonicDeBruijnAdditionMatrix():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Input File
    inFileName = DEBRUIJN_GRAPH_FILENAME
    inFilePath = PATH_TO_DEBRUIJN_GRAPH
    inFile = open(inFilePath + inFileName, 'r')

    outFileName = 'debruijn.germ.addition.matrix'
    outFilePath = PATH_TO_DEBRUIJN_GRAPH
    outFile = open(outFilePath + outFileName, 'w')

    additionMatrix = constructSquareIntegerMatrix(8)

    # Iterate through <inFile>
    for line in inFile:
        # Ignore non-Addition Edges
        lineFields = parseTabSeparatedLine(line)
        if lineFields[2] != 'addition':
            continue

        # Get Germ Layer Hashes
        hashA = getModuleGermLayerHash(lineFields[0])
        hashB = getModuleGermLayerHash(lineFields[1])

        additionMatrix[hashA][hashB] = additionMatrix[hashA][hashB] + 1

    # Write Matrix to File
    writeMatrixToFile(additionMatrix, outFile)

    # Close File and DB Connection
    inFile.close()
    outFile.close()
    cli.close()

    return

# getEmbryonicDeBruijnRemovalMatrix: Creates matrix representing
# the number of removal similarities between vertices of two types. The
# "type" is the germ layer hash of the modules representing the two nodes.
# i.e. (MESO)-->(MESO+ENDO) increments M[2][6] by 1.
def getEmbryonicDeBruijnRemovalMatrix():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Input File
    inFileName = DEBRUIJN_GRAPH_FILENAME
    inFilePath = PATH_TO_DEBRUIJN_GRAPH
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    outFileName = 'debruijn.germ.removal.matrix'
    outFilePath = PATH_TO_DEBRUIJN_GRAPH
    outFile = open(outFilePath + outFileName, 'w')

    removalMatrix = constructSquareIntegerMatrix(8)

    # Iterate through <inFile>
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        # Ignore non-Removal Edges
        if lineFields[2] != 'removal':
            continue

        # Get Germ Layer Hashes
        hashA = getModuleGermLayerHash(lineFields[0])
        hashB = getModuleGermLayerHash(lineFields[1])

        removalMatrix[hashA][hashB] = removalMatrix[hashA][hashB] + 1

    # Write Matrix to File
    writeMatrixToFile(removalMatrix, outFile)

    # Close File and DB Connection
    inFile.close()
    outFile.close()
    cli.close()

    return
    
# createEmbryonicDeBruijnExchangeMatrix: Creates matrix representing
# the number of exchange similarities between vertices of two types. The
# "type" is the germ layer hash of the modules representing the two nodes.
# i.e. (MESO)-->(MESO+ENDO) increments M[2][6] by 1.
def getEmbryonicDeBruijnExchangeMatrix(shuffle = True):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Input File
    if shuffle:
        inFileName = 'shuffle.' + DEBRUIJN_GRAPH_FILENAME
    else:
        inFileName = DEBRUIJN_GRAPH_FILENAME
    inFilePath = PATH_TO_DEBRUIJN_GRAPH
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    if shuffle:
        outFileName = 'shuffle.debruijn.germ.exchange.matrix'
    else: 
        outFileName = 'debruijn.germ.exchange.matrix'
    outFilePath = PATH_TO_DEBRUIJN_GRAPH
    outFile = open(outFilePath + outFileName, 'w')

    exchangeMatrix = constructSquareIntegerMatrix(8)

    # Iterate through <inFile>
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        # Ignore non-Exchange Edges
        if lineFields[2] != 'exchange':
            continue

        # Get Germ Layer Hashes
        if shuffle:
            hashA = getModuleGermLayerHash(lineFields[0], shuffle = True)
            hashB = getModuleGermLayerHash(lineFields[1])
        else:
            hashA = getModuleGermLayerHash(lineFields[0], shuffle = True)
            hashB = getModuleGermLayerHash(lineFields[1])

        exchangeMatrix[hashA][hashB] = exchangeMatrix[hashA][hashB] + 1

    # Write Matrix to File
    writeMatrixToFile(exchangeMatrix, outFile)

    # Close File and DB Connection
    inFile.close()
    outFile.close()
    cli.close()

    return

# getDeBruijnEdgeTypeCounts: Returns # of Exchanges, Additions, and Removals.
def getDeBruijnEdgeTypeCounts(shuffle = False):
    # Open Input File
    if shuffle:
        inFileName = 'shuffle.' + DEBRUIJN_GRAPH_FILENAME
    else:
        inFileName = DEBRUIJN_GRAPH_FILENAME
    inFilePath = PATH_TO_DEBRUIJN_GRAPH
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    if shuffle:
        outFileName = 'shuffle.debruijn.edge.type.counts'
    else:
        outFileName = 'debruijn.edge.type.counts'
    outFilePath = PATH_TO_DEBRUIJN_GRAPH
    outFile = open(outFilePath + outFileName, 'w')

    numAdditions = 0
    numRemovals = 0
    numExchanges = 0

    # Iterate through Graph
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        if lineFields[2] == 'addition':
            numAdditions = numAdditions + 1
        if lineFields[2] == 'removal':
            numRemovals = numRemovals + 1
        if lineFields[2] == 'exchange':
            numExchanges = numExchanges + 1

    # Write Output
    outFile.write('additions: %d\n' % numAdditions)
    outFile.write('removals: %d\n' % numRemovals)
    outFile.write('exchanges: %d\n' % numExchanges)

    # Close Files
    inFile.close()
    outFile.close()

    return

def getSameOrganSystemEdgeCounts():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Input File
    inFileName = DEBRUIJN_GRAPH_FILENAME
    inFilePath = PATH_TO_DEBRUIJN_GRAPH
    inFile = open(inFilePath + inFileName, 'r')

    numSame = 0
    numDiff = 0

    skipLine = True

    # Iterate through Graph
    for line in inFile:
        if skipLine:
            skipLine = False
            continue

        lineFields = parseTabSeparatedLine(line)

        moduleA = lineFields[0]
        moduleB = lineFields[1]

        if moduleA == moduleB:
            continue

        recordA = modDB.find_one( { 'module_id' : moduleA } )
        recordB = modDB.find_one( { 'module_id' : moduleB } )

        tissuesA = recordA.get('tissue_list')
        if len(tissuesA) > 1:
            continue
        tissuesB = recordB.get('tissue_list')
        if len(tissuesB) > 1:
            continue

        if tissuesA[0] == 'global' or tissuesB[0] == 'global':
            continue

        if tissuesA[0] == 'intersection' or tissuesB[0] == 'intersection':
            numSame = numSame + 1
            continue

        organA = TISSUE_SYSTEM_MAP.get(tissuesA[0])
        organB = TISSUE_SYSTEM_MAP.get(tissuesB[0])

        if organA == organB:
            numSame = numSame + 1
        else:
            numDiff = numDiff + 1

    print numSame
    print numDiff

#createDeBruijnSetGraph(shuffleVal = True)
