#! /usr/bin/python

# geneTissueMap.py
# Author: Santhosh Balasubramanian
# Created: January 27, 2013
# Last Modified: March 24, 2013


# Python Imports
import os

# Library Imports
from pymongo import *

# Global Imports
from settings import *


# - - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

GENE_TISSUE_MAP_NOMENCLATURE = 'primary_gene_id'

# - - - - - - - - - - - GENE:TISSUE MAP CREATION - - - - - - - - - - #


# createGeneTissueMap: Creates <geneTissueMap> from <normTissueGeneMap>.
def createGeneTissueMap(minimum = 0, maximum = 10000):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ntgm = db.normTissueGeneMap
    geneTissueMap = db.geneTissueMap

    firstTissue = True

    # Get Gene List
    tissueRecord = ntgm.find_one()
    geneDict = tissueRecord.get(GENE_TISSUE_MAP_NOMENCLATURE)
    geneList = geneDict.keys()[minimum:maximum]

    # Iterate through Genes
    for gene in geneList:
        if PRINT_PROGRESS:
            print gene
            
        record = {}
        record.update( { GENE_TISSUE_MAP_NOMENCLATURE : gene } )

        for tissue in FUNCTIONAL_TISSUE_LIST:
            tissueRecord = ntgm.find_one( { 'tissue' : tissue } )
            expression = tissueRecord.get(GENE_TISSUE_MAP_NOMENCLATURE).get(gene)
            record.update( { tissue : expression } )

        geneTissueMap.save(record)

    # Close DB Connection
    cli.close()

    return


# normalizeGeneTissueMap: Creates <normGeneTissueMap> from <geneTissueMap>.
def normalizeGeneTissueMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap
    geneTissueMap = db.geneTissueMap

    # Iterate through <geneTissueMap>
    for geneRecord in geneTissueMap.find():
        geneID = geneRecord.get(GENE_TISSUE_MAP_NOMENCLATURE)
        if PRINT_PROGRESS:
            print geneID

        record = { 
            'gene_id' : geneID,
            'tissue_list' : [],
            }

        for tissue, expression in geneRecord.iteritems():
            if tissue == '_id':
                continue
            if tissue == GENE_TISSUE_MAP_NOMENCLATURE:
                continue

            if expression > ARRAY_EXPRESSION_THRESHOLD:
                record['tissue_list'].append(tissue)

        ngtmDB.save(record)

    # Close DB Connection
    cli.close()
    
    return
        

# - - - - - - - - - - ANALYSIS - - - - - - - - - - #


# getGeneTissueHistogram: Writes histogram of # of tissues by # genes.
def getGeneTissueHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneTissueMap = db.geneTissueMap

    # Open Output File
    outFile = open(PATH_TO_EXPRESSION_ANALYSIS + 'gene.tissue.histogram', 'w')

    histList = [0]*(len(FUNCTIONAL_TISSUE_LIST) + 1)

    # Iterate through <geneTissueMap>
    for record in geneTissueMap.find():
        count = 0
        # Iterate through Tissues
        for key, value in record.iteritems():
            if key == 'primary_gene_id':
                continue
            if key == '_id':
                continue

            if value > ARRAY_EXPRESSION_THRESHOLD:
                count = count + 1

        for i in range(len(histList)):
            outFile.write('%d\t%d\n' % (i, histList[i]))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return

# getGeneOrganSystemHistogram: Writes histogram of # of systems by # genes.
def getGeneOrganSystemHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneTissueMap = db.geneTissueMap

    # Open Output File
    outFile = open(PATH_TO_EXPRESSION_ANALYSIS + 'gene.organ_system.histogram', 'w')
    
    histList = [0]*(len(FUNCTIONAL_TISSUE_LIST) + 1)

    # Iterate through <geneTissueMap>
    for record in geneTissueMap.find():
        geneOrganDict = {}
        for organSystem in CANONICAL_ORGAN_SYSTEM_LIST:
            geneOrganDict.update( { organSystem : False } )

        # Iterate through Tissues
        for tissue, value in record.iteritems():
            if tissue == 'primary_gene_id':
                continue
            if tissue == '_id':
                continue

            organ = TISSUE_SYSTEM_MAP.get(tissue)

            if value > ARRAY_EXPRESSION_THRESHOLD:
                geneOrganDict[organ] = True

        count = 0
        for organ, boolean in geneOrganDict.iteritems():
            if boolean:
                count = count + 1

        histList[count] = histList[count] + 1

    for i in range(len(histList)):
        outFile.write('%d\t%d\n' % (i, histList[i]))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
