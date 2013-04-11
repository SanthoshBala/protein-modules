#! /usr/bin/python

# deBruijnSequences.py
# Author: Santhosh Balasubramanian
# Created: March 31, 2013
# Last Modified: March 31, 2013


# Python Imports
from pymongo import *

# Global Imports
from settings import *

# Nomenclature Imports
from nomenclature.proteinID import *

# Utility Imports
from common.strings import *


# getExchangeSimilarGenes: Creates file containing list of genes that
# replace one another in a module, as well as file with pairwise relations.
def getExchangeSimilarGenes():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Input File
    inFileName = 'debruijn.similarity.graph'
    inFilePath = PATH_TO_MODULE_SIMILARITY + inFileName
    inFile = open(inFilePath, 'r')

    # Open Output File
    genePairFileName = 'debruijn.exchange.gene.pairs'
    geneListFileName = 'debruijn.exchange.genes'
    outFilePath = PATH_TO_MODULE_SIMILARITY
    genePairFile = open(outFilePath + genePairFileName, 'w')
    geneListFile = open(outFilePath + geneListFileName, 'w')

    geneSet = set()

    # Collect Gene Set
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)

        # Ignore non-Exchange Relationships
        if lineFields[2] != 'substitution':
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
    genePairFileName = 'debruijn.exchange.gene.pairs'
    inFilePath = PATH_TO_MODULE_SIMILARITY
    genePairFile = open(inFilePath + genePairFileName, 'r')

    # Open Output Files
    proteinListName = 'debruijn.exchange.proteins'
    proteinPairName = 'debruijn.exchange.protein.pairs'
    outFilePath = PATH_TO_MODULE_SIMILARITY
    proteinListFile = open(outFilePath + proteinListName, 'w')
    proteinPairFile = open(outFilePath + proteinPairName, 'w')

    geneSet = set()
    proteinSet = set()

    # Collect Protein Set
    for line in genePairFile:
        lineFields = parseTabSeparatedLine(line)

        geneA = lineFields[0]
        geneB = lineFields[1]

        # Get Proteins
        if 'ENSG' in geneA:
            proteinA = getRefseqProteinForEnsemblGene(geneA)
        else:
            proteinA = getRefseqProteinForEntrezGene(geneA)
        
        if not proteinA:
            continue

        if 'ENSG' in geneB:
            proteinB = getRefseqProteinForEnsemblGene(geneB)
        else:
            proteinB = getRefseqProteinForEntrezGene(geneB)
            
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
