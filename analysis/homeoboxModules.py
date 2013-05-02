#! /usr/bin/python

# homeoboxModules.py
# Author: Santhosh Balasubramanian
# Created: April 17, 2013
# Last Modified: April 28, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphDraw import *

# Nomenclature Imports
from nomenclature.geneID import *

# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

HOMEOBOX_MODULE_FILENAME = 'homeobox.module.list'
PATH_TO_HOMEOBOX_MODULES = PATH_TO_ANALYSIS + 'homeobox/'

HOMEOBOX_MODULE_IDS = [
    'MOD960600003255',
    'MOD960600004406',
    'MOD960600001185',
    'MOD960600001325',
    'MOD960600005320',
    'MOD960600006998',
    'MOD960600004578',
    'MOD960600006433',
    'MOD960600006344',
    'MOD960600001954'
]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# drawHomeoboxModules
def drawHomeoboxModules():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Iterate through moduleIDs
    for moduleID in HOMEOBOX_MODULE_IDS:
        if PRINT_PROGRESS:
            print moduleID

        # Draw Module using Gene Symbol
        drawModule(moduleID, symbol = True)

        # Close DB Connection
    cli.close()

    return

# writeHomeoboxModules
def writeHomeoboxModules():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFileName = HOMEOBOX_MODULE_FILENAME
    outFilePath = PATH_TO_HOMEOBOX_MODULES
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through Modules
    for moduleID in HOMEOBOX_MODULE_IDS:
        if PRINT_PROGRESS:
            print moduleID

        outFile.write('%s\n' % moduleID)

        record = modDB.find_one( { 'module_id' : moduleID } )

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
