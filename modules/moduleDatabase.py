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


# - - - - - - - - - - RAW MODULE FILE CREATION - - - - - - - - - - #


# Raw module files just contain a list of genes in each row, as is the
# raw output of SPICi. These are later converted to files of just module_ids.

# createAllRawModuleFiles: Creates raw module files for all tissues.
def createAllRawModuleFiles():
    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue
        createTissueRawModuleFile(tissue)

    return

# createTissueRawModuleFile: Creates raw module file for <tissue>.
def createTissueRawModuleFile(tissue):
    # Get Filename for Tissue Graph
    inFilename = GENEMANIA_TISSUE_BASE_FILENAME % \
        (tissue, DESIRED_INTERACTION_TYPE)

    # Convert Genemania to SPICi Format
    canonicalizeGenemaniaGraph(inFilename)
    canonFilename = CANON_GENEMANIA_TISSUE_BASE_FILENAME % \
        (tissue, DESIRED_INTERACTION_TYPE)
    inFilePath = PATH_TO_TISSUE_SUBGRAPHS + canonFilename
    
    # Open Raw Module File Path
    outFilename = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
    outFilePath = PATH_TO_TISSUE_SPICI_MODULES + outFilename

    # Run SPICi
    command = [ PATH_TO_SPICI_BINARY, '-i', inFilePath, '-o', outFilePath, 
                '-s', str(MODULE_SIZE_THRESHOLD),'-m', '2']
    subprocess.call(command)

    # Remove canonical input file
    os.chdir(PATH_TO_TISSUE_SUBGRAPHS)
    os.remove(canonFilename)

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
def createModuleDB:
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Create Number System for Modules: MOD960600000000
    moduleNum = 0

    # Fill DB with Modules from Tissue Subgraphs
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue

        # Open Tissue Modules File
        inFileName = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
        inFilePath = PATH_TO_TISSUE_SPICI_MODULES + inFileName
        inFile = open(inFilePath, 'r')

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
def createModuleGermLayerField():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Iterate through DB
    for moduleRecord in modDB.find(snapshot = True):
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


# createModuleGeneUniversalityField: Adds 'gene_universality' to modDB. 
# 'gene_universality' is float representing % of tissues holding genes.
def createModuleGeneUniversalityField():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    ngtmDB = db.normGeneTissueMap

    # Iterate through DB
    for moduleRecord in modDB.find(snapshot = True):
        moduleID = module.get('module_id')
        
        if PRINT_PROGRESS:
            print moduleID

        geneList = module.get('gene_list')
        numGenes = len(geneList)

        # Get Mean # of Tissues
        tissueSum = 0
        for geneID in geneList:
            geneRecord = ngtmDB.find_one( { 'primary_gene_id' : geneID } )
            if not geneRecord:
                continue
            numTissues = len(geneRecord.get('tissue_list'))
            tissueSum = tissueSum + numTissues
        tissuesPerGene = tissueSum/float(numGenes)
        geneUniversality = tissuesPerGene/float(len(FUNCTIONAL_TISSUE_LIST))
        
        # Update <modDB>
        moduleRecord.update( { 'gene_universality' : geneUniversality } )
        modDB.save(moduleRecord)

    # Close DB Connection
    cli.close()
    
    return


# createModuleHousekeepingIndexField: Adds 'housekeeping' to modDB. 'housekeeping'
# represents the % of genes in module that are housekeeping genes.
def createModuleHousekeepingIndexField():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    ngtmDB = db.normGeneTissueMap

    # Iterate through DB
    for moduleRecord in modDB.find(snapshot = True):
        moduleID = module.get('module_id')
        if PRINT_PROGRESS:
            print moduleID

        geneList = module.get('gene_list')
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
def constructGeneModuleMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    geneModMap = db.geneModuleMap

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
def createAllModuleIDFiles():
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue
        createTissueModuleIDFile(tissue)

# createTissueModuleIDFile: Creates module ID file for <tissue>.
def createTissueModuleIDFile(tissue):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Raw File
    inFileName = TISSUE_SPICI_MODULES_BASE_FILENAME % tissue
    inFile = open(PATH_TO_SPICI_MODULES + inFileName, 'r')

    # Open Output File
    outFileName = TISSUE_MODULE_IDS_BASE_FILENAME % tissue
    outFile = open(PATH_TO_TISSUE_MODULE_IDS + outFileName, 'w')

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
    
