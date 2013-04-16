#! /usr/bin/python

# moduleDatabase.py
# Author: Santhosh Balasubramanian
# Created: February 16, 2013
# Last Modified: March 26, 2013


# Python Imports
import os
import subprocess

# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graphs Imports
from graphs.graphIO import *
from graphs.graphUtil import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


SHUFFLE_SPICI_MODULES_BASE_FILENAME = 'shuffle.%s.modules'
PATH_TO_SHUFFLE_SPICI_MODULES = PATH_TO_MODULES + 'shuffle_spici_modules/'
SHUFFLE_MODULE_IDS_BASE_FILENAME = 'shuffle.%s.module.ids'
PATH_TO_SHUFFLE_MODULE_IDS = PATH_TO_MODULES + 'shuffle_module_ids/'


# - - - - - - - - - - RAW MODULE FILE CREATION - - - - - - - - - - #


# Raw module files just contain a list of genes in each row, as is the
# raw output of SPICi. These are later converted to files of just module_ids.

# createAllRawModuleFiles: Creates raw module files for all tissues.
def createAllRawModuleFiles(minimum = 0, maximum = 85, shuffle = False):
    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue
        createTissueRawModuleFile(tissue, shuffle = shuffle)

    return

# createTissueRawModuleFile: Creates raw module file for <tissue>.
def createTissueRawModuleFile(tissue, shuffle = False):
    # Get Filename for Tissue Graph
    if shuffle:
        inFilePath = PATH_TO_SHUFFLE_TISSUE_SUBGRAPHS
        inFileName = SHUFFLE_TISSUE_SUBGRAPH_BASE_FILENAME % tissue
    else:
        inFilePath = PATH_TO_TISSUE_SUBGRAPHS
        inFileName = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, 'ppi')

    # Convert Genemania to SPICi Format
    canonicalizeGenemaniaGraph(inFilePath, inFileName)
    canonFileName = 'canon.' + inFileName
    inFile = inFilePath + canonFileName
    
    # Open Raw Module File Path
    if shuffle:
        outFileName = SHUFFLE_SPICI_MODULES_BASE_FILENAME % tissue
        outFilePath = PATH_TO_SHUFFLE_SPICI_MODULES + outFileName
    else:
        outFileName = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
        outFilePath = PATH_TO_TISSUE_SPICI_MODULES + outFileName

    # Run SPICi
    command = [ PATH_TO_SPICI_BINARY, '-i', inFile, '-o', outFilePath, 
                '-s', str(MODULE_SIZE_THRESHOLD),'-m', '2']
    subprocess.call(command)

    # Remove canonical input file
    os.chdir(inFilePath)
    os.remove(canonFileName)

    return

# createGraphRawModuleFile: Creates raw module file for <graph>.
def createGraphRawModuleFile(graph, outputPath):
    graphTempFile = PATH_TO_NETWORKS + 'graph.temp'

    # Write Graph to Disk
    writeCanonicalGraphToDisk(graph, graphTempFile)

    # Run SPICi
    command = [ PATH_TO_SPICI_BINARY, '-i', graphTempFile, '-o', outputPath, 
                '-s', str(MODULE_SIZE_THRESHOLD), '-m', '2']
    subprocess.call(command)
    
    # Remove Graph File
    os.remove(graphTempFile)

    return


# - - - - - - - - - - MODULE DB CREATION - - - - - - - - - - #


# createModuleDB: Creates <modules> a database for all modules.
def createModuleDB(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
    else:
        modDB = db.modules

    # Create Number System for Modules: MOD960600000000
    moduleNum = 0

    # Fill DB with Modules from Tissue Subgraphs
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue

        # Open Tissue Modules File
        if shuffle:
            inFileName = SHUFFLE_SPICI_MODULES_BASE_FILENAME % tissue
            inFilePath = PATH_TO_SHUFFLE_SPICI_MODULES
        else:
            inFileName = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
            inFilePath = PATH_TO_TISSUE_SPICI_MODULES
        inFile = open(inFilePath + inFileName, 'r')

        # Iterate through File
        for line in inFile:
            geneList = parseTabSeparatedLine(line)
            # Sort to get Canonical Representation
            geneList.sort()

            # Try to Find Module
            record = modDB.find_one( { 'gene_list' : geneList } )
            if record:
                if tissue in record['tissue_list']:
                    continue
                else:
                    record['tissue_list'].append(tissue)
                    modDB.save(record)
            else:
                # Create New Module ID
                suffixNum = 960600000000 + moduleNum
                moduleNum = moduleNum + 1
                moduleID = 'MOD%d' % suffixNum

                # Add to <modDB>
                modDB.insert( {
                        'gene_list' : geneList,
                        'tissue_list' : [ tissue ],
                        'module_id' : moduleID
                        } )

        # Close File
        inFile.close()
        
    # Close DB Connection
    cli.close()

    return

# createModuleGermLayerField: Adds 'germ_layer' field to modDB.
def createModuleGermLayerField(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.modules
    else:
        modDB = db.shuffleModules

    # Iterate through DB
    for moduleRecord in modDB.find(snapshot = True, timeout = False):
        # Get Germ Layers for <moduleRecord> Tissues
        tissueList = moduleRecord.get('tissue_list')
        germLayerSet = set()
        for tissue in tissueList:
            if tissue == 'global' or tissue == 'intersection':
                continue
            germLayer = TISSUE_GERM_LAYER_MAP.get(tissue)
            germLayerSet.add(germLayer)
            
        moduleRecord.update( { 'germ_layer' : list( germLayerSet ) } )
        
        modDB.save(moduleRecord)
        
    # Close DB Connection
    cli.close()

    return


# createModuleProteinUniversalityField: Adds 'protein_universality' to modDB. 
# 'protein_universality' is float representing % of tissues holding genes.
def createModuleProteinUniversalityField(skipVal = 0, limitVal = 0, 
                                         shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
        gtmDB = db.shuffleNormGeneTissueMap
    else:
        modDB = db.modules
        gtmDB = db.normGeneTissueMap

    # Iterate through DB
    for moduleRecord in modDB.find(skip = skipVal, limit = limitVal, 
                                   snapshot = True, timeout = False):
        moduleID = moduleRecord.get('module_id')
        
        if PRINT_PROGRESS:
            print moduleID

        geneList = moduleRecord.get('gene_list')
        numGenes = len(geneList)

        # Get Mean # of Tissues
        tissueSum = 0
        for geneID in geneList:
            geneRecord = gtmDB.find_one( { 'primary_gene_id' : geneID } )
            if not geneRecord:
                continue
            numTissues = len(geneRecord.get('tissue_list'))
            tissueSum = tissueSum + numTissues
        tissuesPerGene = tissueSum/float(numGenes)
        proteinUniversality = tissuesPerGene/float(len(FUNCTIONAL_TISSUE_LIST))
        
        # Update <modDB>
        moduleRecord.update( { 'protein_universality' : proteinUniversality } )
        modDB.save(moduleRecord)

    # Close DB Connection
    cli.close()
    
    return


# createModuleHousekeepingIndexField: Adds 'housekeeping' to modDB. 
# 'housekeeping' represents the % of genes in module that are housekeeping.
def createModuleHousekeepingIndexField(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
        gtmDB = db.shuffleGeneTissueMap
    else:
        modDB = db.modules
        gtmDB = db.normGeneTissueMap

    # Iterate through DB
    for moduleRecord in modDB.find(snapshot = True, timeout = False):
        moduleID = moduleRecord.get('module_id')
        if PRINT_PROGRESS:
            print moduleID

        geneList = moduleRecord.get('gene_list')
        numGenes = len(geneList)
        numHousekeeping = 0

        if numGenes < MODULE_SIZE_THRESHOLD:
            continue

        # Iterate through Genes
        for geneID in geneList:
            geneRecord = ngtmDB.find_one( { 'primary_gene_id' : geneID } )
            if geneRecord:
               tissueList = geneRecord.get('tissue_list')
               numTissues = len(tissueList)
               if numTissues == len(FUNCTIONAL_TISSUE_LIST):
                   numHousekeeping = numHousekeeping + 1

        # Update DB
        housekeepingIndex = numHousekeeping/float(numGenes)
        moduleRecord.update( { 'housekeeping' : housekeepingIndex } )
        modDB.save(moduleRecord)

    # Close DB Connection
    cli.close()

    return

# createGeneModuleMap: Creates <geneModuleMap>, which contains modules to
# which each gene belongs in each tissue.
def constructGeneModuleMap(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.modules
        geneModMap = db.geneModuleMap
    else:
        modDB = db.shuffleModules
        geneModMap = db.shuffleGeneModuleMap

    # Iterate through Modules
    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        tissueList = moduleRecord.get('tissue_list')
        geneList = moduleRecord.get('gene_list')

        # Iterate through <geneList>
        for geneID in geneList:
            # Get or Create Record for <gene>
            geneRecord = geneModMap.find_one( { 'gene_id' : geneID } )
            if not geneRecord:
                geneRecord = { 'gene_id' : geneID }

            # Add Tissues to Record
            for tissue in tissueList:
                geneRecord.update( { tissue : moduleID } )

            geneModMap.save(geneRecord)

    # Close DB Connection
    cli.close()

    return


# - - - - - - - - - - MODULE ID FILE CREATION - - - - - - - - - - #


# Module ID files contain a list of module IDs, with one per line.

# createAllModuleIDFiles: Creates module ID files for all tissues.
def createAllModuleIDFiles(shuffleVal = False):
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue
        createTissueModuleIDFile(tissue, shuffle = shuffleVal)

# createTissueModuleIDFile: Creates module ID file for <tissue>.
def createTissueModuleIDFile(tissue, shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        modDB = db.shuffleModules
    else:
        modDB = db.modules

    # Open Raw File
    if shuffle:
        inFilePath = PATH_TO_SHUFFLE_MODULES
        inFileName = SHUFFLE_SPICI_MODULES_BASE_FILENAME % tissue
    else:
        inFilePath = PATH_TO_TISSUE_SPICI_MODULES
        inFileName = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    if shuffle:
        outFilePath = PATH_TO_SHUFFLE_MODULE_IDS
        outFileName = SHUFFLE_MODULE_IDS_BASE_FILENAME % tissue
    else:
        outFilePath = PATH_TO_TISSUE_MODULE_IDS
        outFileName = TISSUE_MODULE_IDS_BASE_FILENAME % tissue
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through Raw File
    for line in inFile:
        geneList = parseTabSeparatedLine(line)
        # Sort to get Canonical Representation
        geneList.sort()

        # Get Record for Module
        moduleRecord = modDB.find_one( { 'gene_list' : geneList } )
        moduleID = moduleRecord.get('module_id')

        # Write Module ID
        outFile.write('%s\n' % moduleID)

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
    
