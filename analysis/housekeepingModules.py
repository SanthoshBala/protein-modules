#! /usr/bin/python

# housekeepingModules.py
# Author: Santhosh Balasubramanian
# Created: April 17, 2013
# Last Modified: April 17, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphDraw import *

# Nomenclature Imports
from nomenclature.geneID import *

# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

HOUSEKEEPING_MODULE_FILENAME = 'housekeeping.module.list'
PATH_TO_HOUSEKEEPING_MODULES = PATH_TO_ANALYSIS + 'housekeeping/'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# drawHousekeepingModules
def drawHousekeepingModules():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    moduleRecords = modDB.find( { 'tissue_list' : 'intersection' } )
    
    # Iterate through <moduleRecords>
    for record in moduleRecords:
        moduleID = record.get('module_id')
        if PRINT_PROGRESS:
            print moduleID

        # Draw Module using Gene Symbol
        drawModule(moduleID, symbol = True)

        # Close DB Connection
    cli.close()

    return

# writeHousekeepingModules
def writeHousekeepingModules():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFileName = HOUSEKEEPING_MODULE_FILENAME
    outFilePath = PATH_TO_HOUSEKEEPING_MODULES
    outFile = open(outFilePath + outFileName, 'w')

    moduleRecords = modDB.find( { 'tissue_list' : 'intersection' } )
    
    # Iterate through Modules
    for record in moduleRecords:
        moduleID = record.get('module_id')
        if PRINT_PROGRESS:
            print moduleID

        outFile.write('%s\n' % moduleID)

        # Iterate through Genes
        for geneID in record.get('gene_list'):
            if PRINT_PROGRESS:
                print geneID

            if 'ENSG' in geneID:
                entrezID = getEntrezForEnsemblGene(geneID)
                if not entrezID:
                    symbol = geneID
                else:
                    symbol = getNCBIOfficialSymbol(entrezID)
            else:
                entrezID = geneID
                symbol = getNCBIOfficialSymbol(entrezID)
            # Write Gene Symbols
            outFile.write('\t%s & %s & %s \\\\\n' % (geneID, entrezID, symbol))

        outFile.write('\n')

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
