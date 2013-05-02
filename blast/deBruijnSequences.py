#! /usr/bin/python

# deBruijnSequences.py
# Author: Santhosh Balasubramanian
# Created: March 31, 2013
# Last Modified: March 31, 2013


# Python Imports
from math import *
from random import *

# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Nomenclature Imports
from nomenclature.geneID import *
from nomenclature.proteinID import *

# Utility Imports
from common.strings import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

DEBRUIJN_GRAPH_FILENAME = 'debruijn.similarity.graph'

PATH_TO_EXCHANGE_SIMILARITY = PATH_TO_MODULE_SIMILARITY + 'exchange/'
SIMILAR_GENE_LIST_FILENAME = 'debruijn.exchange.genes'
SIMILAR_GENE_PAIRS_FILENAME = 'debruijn.exchange.gene.pairs'
SIMILAR_PROTEIN_LIST_FILENAME = 'debruijn.exchange.proteins'
SIMILAR_PROTEIN_PAIRS_FILENAME = 'debruijn.exchange.protein.pairs' 
SIMILAR_PROTEIN_SHARD_BASE_FILENAME = 'debruijn.exchange.proteins.shard.%d'
PROTEINS_PER_SHARD = 100


# - - - - - - - - - - FUNCTIONS - - - - - - - - - - #


# getExchangeSimilarGenes: Creates file containing list of genes that
# replace one another in a module, as well as file with pairwise relations.
def getExchangeSimilarGenes(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
    else:
        modDB = db.modules

    # Open Input File
    inFileName = DEBRUIJN_GRAPH_FILENAME
    inFilePath = PATH_TO_MODULE_SIMILARITY
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    genePairFileName = SIMILAR_GENE_PAIRS_FILENAME
    geneListFileName = SIMILAR_GENE_LIST_FILENAME
    outFilePath = PATH_TO_MODULE_SIMILARITY
    genePairFile = open(outFilePath + genePairFileName, 'w')
    geneListFile = open(outFilePath + geneListFileName, 'w')

    geneSet = set()

    # Collect Gene Set
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)

        # Ignore non-Exchange Relationships
        if lineFields[2] != 'exchange':
            continue

        moduleAID = lineFields[0]
        moduleBID = lineFields[1]

        # Get Gene Sets
        moduleARecord = modDB.find_one( { 'module_id' : moduleAID } )
        moduleAGenes = set(moduleARecord.get( 'gene_list' ))
        moduleBRecord = modDB.find_one( { 'module_id' : moduleBID } )
        moduleBGenes = set(moduleBRecord.get( 'gene_list' ))
        exchangeGenes = list(moduleAGenes.symmetric_difference(moduleBGenes))
        
        # Write Pairs
        if exchangeGenes[0] < exchangeGenes[1]:
            geneSet.add(exchangeGenes[0])
            geneSet.add(exchangeGenes[1])
            genePairFile.write('%s\t%s\texchange\n' % (exchangeGenes[0], 
                                                       exchangeGenes[1]))

    # Write List
    for gene in geneSet:
        geneListFile.write('%s\n' % gene)

    # Close Files and DB Connection
    genePairFile.close()
    geneListFile.close()
    cli.close()

    return
    

# getExchangeSimilarProteins: Creates file containing list of proteins that
# replace one another in a module, as well as file with pairwise relations.
def getExchangeSimilarProteins():
    # Open Input Files
    genePairFileName = SIMILAR_GENE_PAIRS_FILENAME
    inFilePath = PATH_TO_EXCHANGE_SIMILARITY
    genePairFile = open(inFilePath + genePairFileName, 'r')

    # Open Output Files
    proteinListFileName = SIMILAR_PROTEIN_LIST_FILENAME
    proteinPairFileName = SIMILAR_PROTEIN_PAIRS_FILENAME
    outFilePath = PATH_TO_EXCHANGE_SIMILARITY
    proteinListFile = open(outFilePath + proteinListFileName, 'w')
    proteinPairFile = open(outFilePath + proteinPairFileName, 'w')

    geneSet = set()
    proteinSet = set()

    # Collect Protein Set
    for line in genePairFile:
        lineFields = parseTabSeparatedLine(line)

        geneA = lineFields[0]
        geneB = lineFields[1]

        # Get Proteins
        if 'ENSG' in geneA:
            proteinA = getUniprotForEnsemblGene(geneA)
        else:
            proteinA = getUniprotForEntrezGene(geneA)
        
        if not proteinA:
            continue

        if 'ENSG' in geneB:
            proteinB = getUniprotForEnsemblGene(geneB)
        else:
            proteinB = getUniprotForEntrezGene(geneB)
            
        if not proteinB:
            continue

        proteinSet.add(proteinA)
        proteinSet.add(proteinB)

        # Write Pairs
        proteinPairFile.write('%s\t%s\texchange\n' % (proteinA, proteinB))

    
    # Write Output
    for protein in proteinSet:
        proteinListFile.write('%s\n' % protein)

    proteinPairFile.close()
    proteinListFile.close()

    return


# shardExchangeSimilarProteinList: BLAST can only accept a limited number
# of queries, so shard the list into N files such that each file only has
# 100 proteins, easing load on NCBI servers
def shardExchangeSimilarProteinList():
    # Open Input File
    inFilePath = PATH_TO_EXCHANGE_SIMILARITY
    inFileName = SIMILAR_PROTEIN_LIST_FILENAME
    inFile = open(inFilePath + inFileName, 'r')

    proteinList = []
    shardNum = 0

    # Iterate through <inFile>
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        protein = lineFields[0]
        proteinList.append(protein)

    numProteins = len(proteinList)

    # Randomize Protein List
    augmentedProteinList = []
    shuffle(proteinList)
    augmentedProteinList.extend(proteinList)
    shuffle(proteinList)
    augmentedProteinList.extend(proteinList)
    numShards = int(ceil(2*numProteins/float(PROTEINS_PER_SHARD)))

    # Shard <proteinList> into Files of Size N
    for i in range(numShards):
        outFileName = SIMILAR_PROTEIN_SHARD_BASE_FILENAME % (i + 1)
        outFilePath = PATH_TO_EXCHANGE_SIMILARITY
        outFile = open(outFilePath + outFileName, 'w')

        proteinSubList = proteinList[i*PROTEINS_PER_SHARD:
                                         (i + 1)*PROTEINS_PER_SHARD]

        for protein in proteinSubList:
            outFile.write('%s\n' % protein)

        outFile.close()

    # Close Input File
    inFile.close()

    return
