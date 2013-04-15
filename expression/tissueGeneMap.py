#! /usr/bin/python

# tissueGeneMap.py
# Author: Santhosh Balasubramanian
# Created: January 17, 2013
# Last Modified: March 24, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.databases import *
from common.statistics import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - - #


TISSUE_GENE_MAP_NOMENCLATURES = [
    'primary_gene_id',
    'ensembl_gene',
    'entrez_gene',
    'gene_symbol',
]


# - - - - - - - - - - TISSUE:GENE MAP CREATION - - - - - - - - - - - #


# createTissueGeneMap: Create <tissueGeneMap> from <tissueProbeMap> and <probeGeneMap>.
def createTissueGeneMap(skipVal = 0, limitVal = 80):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    tissueGeneMap = db.tissueGeneMap
    extProbeGeneMap = db.extProbeGeneMap
    normTissueProbeMap = db.normTissueProbeMap
    
    # Iterate through Tissues
    for tissueRecord in normTissueProbeMap.find(timeout = False, skip = skipVal,
                                          limit = limitVal):
        if PRINT_PROGRESS:
            print tissueRecord.get('tissue')

        # Initialize Record
        record = {
            'tissue' : tissueRecord.get('tissue'),
            'geo_id' : tissueRecord.get('geo_id')
            }
        for nomenclature in TISSUE_GENE_MAP_NOMENCLATURES:
            record.update( { nomenclature : {} } )

        # Iterate through Probes
        for probeID, expression in tissueRecord.get('probes').iteritems():
            # Construct Probe Query
            query = constructQueryForProbe(probeID)
            probeRecord = extProbeGeneMap.find_one(query)
            
            # Iterate through "entrez_id", "unigene_id", etc.
            for nomenclature in TISSUE_GENE_MAP_NOMENCLATURES:
                # i.e. Iterate through all ENSG for probe['ensembl_gene']
                if probeRecord.get(nomenclature) == None:
                    continue
                for gene in probeRecord.get(nomenclature):
                    # If Nomenclature == 'gene_title', get rid of '.'
                    if '.' in gene:
                        gene = ''.join( gene.split('.') )
                        
                    if record[nomenclature].get(gene):
                        record[nomenclature][gene].append(expression)
                    else:
                        record[nomenclature].update( { gene : [ expression ] } )
        tissueGeneMap.save(record)
        
    # Close DB Connection
    cli.close()

    return

# normalizeTissueGeneMap: Normalize <tissueGeneMap> by taking median of probe values.
def normalizeTissueGeneMap(skipVal = 0, limitVal = 80):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    tissueGeneMap = db.tissueGeneMap
    normTissueGeneMap = db.normTissueGeneMap

    # Iterate through Tissues
    for tissueRecord in tissueGeneMap.find(skip=skipVal, limit=limitVal):
        record = {}
        record.update( { 'tissue' : tissueRecord.get('tissue') } )
        record.update( { 'geo_id' : tissueRecord.get('geo_id') } )

        # Iterate through Genes and Normalize
        for nomenclature in TISSUE_GENE_MAP_NOMENCLATURES:
            record.update( { nomenclature : {} } )
            nomenclatureDict = tissueRecord.get(nomenclature)

            for gene, expressionVals in nomenclatureDict.iteritems():
                value = median(expressionVals)
                record[nomenclature].update( { gene : value } )
                
        normTissueGeneMap.save(record)

    # Close DB Connection
    cli.close()
    
    return


# shuffleNormTissueGeneMap: Uses shuffleGeneTissueMap to reverse engineer
# shuffleNormTissueGeneMap.
def shuffleNormTissueGeneMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ntgmDB = db.shuffleNormTissueGeneMap
    shuffleGTM = db.shuffleGeneTissueMap

    # Find All Genes for Tissue Above Threshold Expression
    for tissue in FUNCTIONAL_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue
            
        expressedGenes = shuffleGTM.find( { tissue :
                                                { '$gt' : 
                                                  ARRAY_EXPRESSION_THRESHOLD }
                                            } )

        tissueRecord = {
            'tissue' : tissue,
            'primary_gene_id' : {}
            }
        
        nomenclatureDict = tissueRecord.get('primary_gene_id')
        
        for geneRecord in expressedGenes:
            geneID = geneRecord.get('primary_gene_id')
            nomenclatureDict.update( { geneID : ARRAY_EXPRESSION_THRESHOLD + 1} )

        ntgmDB.save(tissueRecord)

    # Close DB Connection
    cli.close()

    return


# - - - - - - - - - - ANALYSIS - - - - - - - - - - - #


# analyzeTissueGeneMap: Finds # of expressed, unexpressed, and ambiguous genes
# in each tissue, by nomenclature.
def analyzeTissueGeneMap():
    # Open DB connection
    cli = MongoClient()
    db = cli.db
    tissueGeneMap = db.tissueGeneMap

    # Open Output File
    outFile = open(PATH_TO_EXPRESSION_ANALYSIS + 'tissue.gene.map.analysis', 'w')
    outFile.write('Analyzing Tissue:Gene Map...\n\n')

    # Iterate through Tissues
    for tissueRecord in tissueGeneMap.find():
        outFile.write('Tissue: %s\n' % tissueRecord.get('tissue'))
        outFile.write('\tSamples: %s\n' % str(tissueRecord.get('geo_ids')))
        outFile.write('\tNomenclature\tExpressed\tUnexpressed\tAmbiguous\n')
        
        # Iterate through Nomenclatures
        for nomenclature in TISSUE_GENE_MAP_NOMENCLATURES:
            expressedGenes = set()
            unexpressedGenes = set()
            ambiguousGenes = set()
            nomenclatureDict = tissueRecord.get(nomenclature)

            # Iterate through Genes
            for gene, expressionList in nomenclatureDict.iteritems():
                numExpressed = 0
                for expression in expressionList:
                    if expression > ARRAY_EXPRESSION_THRESHOLD:
                        numExpressed = numExpressed + 1
                        
                if numExpressed == len(expressionList):
                    expressedGenes.add(gene)
                elif numExpressed == 0:
                    unexpressedGenes.add(gene)
                else:
                    ambiguousGenes.add(gene)
        
            expCount = len(expressedGenes)
            unexpCount = len(unexpressedGenes)
            ambCount = len(ambiguousGenes)
            outFile.write('\t%s\t%s\t%s\t%s\n' % \
                        (nomenclature, expCount, unexpCount, ambCount))
        outFile.write('\n')

    # Close File and DB Connection
    outFile.close()
    cli.close()
    
    return
