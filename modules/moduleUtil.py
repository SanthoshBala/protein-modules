#! /usr/bin/python

# moduleUtil.py
# Author: Santhosh Balasubramanian
# Created: February 16, 2013
# Last Modified: March 25, 2013


# Python Imports
import subprocess
import os

# Global Imports
from settings import *

# Utility Imports
from common.strings import *

# Graphs Imports
from graphs.graphIO import *
from graphs.graphUtil import *


# - - - - - - - - - - MODULE RETRIEVAL - - - - - - - - - - #


# getTissueModuleIDSet: Returns IDs for all modules found in <tissue>.
def getTissueModuleIDSet(tissue):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Get Record List
    recordList = modDB.find( { 'tissue_list' : tissue } )

    # Get Module IDs
    moduleIDList = [ record.get('module_id') for record in recordList ]

    return set(moduleIDList)

# getTissueRawModuleSet: Returns modules found in <tissue>, as defined by
# the raw gene lists generated by SPICi.
def getTissueRawModuleSet(tissue):
    # Open Raw Module File
    inFileName = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
    inFile = open(PATH_TO_TISSUE_SPICI_MODULES + inFileName, 'r')

    # Initialize <moduleSet>
    moduleSet = set()
    
    # Iterate through File
    for line in inFile:
        geneList = parseTabSeparatedLine(line)
        geneList.sort()

        geneTuple = tuple(geneList)

        moduleSet.add(geneTuple)
    
    # Close File
    inFile.close()

    return moduleSet

# getModuleID: Returns module ID corresponding to <geneList>.
def getModuleID(geneList):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Verify <geneList> is Sorted
    geneList.sort()
    record = modDB.find_one( { 'gene_list' : geneList } )
    
    if record:
        return record.get('module_id')
    else:
        return None

# getTissueModularizedGeneSet: Returns all genes in modules in <tissue>.
def getTissueModularizedGeneSet(tissue):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Get Record List
    recordList = modDB.find( { 'tissue_list' : tissue } )

    geneSet = set()

    # Get Genes
    for record in recordList:
        geneList = record.get('gene_list')
        for gene in geneList:
            geneSet.add(gene)

    # Return Gene Set
    return geneSet


# - - - - - - - - - - MODULE PROPERTIES - - - - - - - - - - #


# getModuleGermLayerHash: Gets germ layer hash for <moduleID>. The germ layer
# hash describes, for a given module, which germ layers it is found in.
def getModuleGermLayerHash(moduleID):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Get Germ Layers
    moduleRecord = modDB.find_one( { 'module_id' : moduleID } )
    germLayers = moduleRecord.get('germ_layer')

    value = 0

    # Compute Hash
    if len(germLayers) == 1:
        if 'endoderm' in germLayers:
            value = 1
        if 'mesoderm' in germLayers:
            value = 2
        if 'ectoderm' in germLayers:
            value = 3
    if len(germLayers) == 2:
        if 'endoderm' not in germLayers:
            value = 6
        if 'mesoderm' not in germLayers:
            value = 5
        if 'ectoderm' not in germLayers:
            value = 4
    if len(germLayers) == 3:
        value = 7
        
    return value
