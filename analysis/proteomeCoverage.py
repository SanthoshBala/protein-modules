#! /usr/bin/python

# proteomeCoverage.py
# Author: Santhosh Balasubramanian
# Created: April 14, 2013
# Last Modified: April 14, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphIO import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

PATH_TO_PROTEOME_COVERAGE = PATH_TO_ANALYSIS + 'proteome_coverage/'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# getProteomeCoverage: Creates file containing number of proteins per tissue
# according to GeneMANIA and Su et al. data.
def getProteomeCoverage(shuffle = False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    if shuffle:
        ngtmDB = db.shuffleNormGeneTissueMap
    else:
        ngtmDB = db.normGeneTissueMap

    # Open Output File
    if shuffle:
        outFileName = 'shuffle.proteome.coverage'
        outFilePath = PATH_TO_PROTEOME_COVERAGE
    else:
        outFileName = 'proteome.coverage'
        outFilePath = PATH_TO_PROTEOME_COVERAGE
    outFile = open(outFilePath + outFileName, 'w')

    outFile.write('Tissue\tMicroarray\tPPI\nCoverage')

    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue

        # Get Num Proteins According to Microarray
        if tissue == 'global':
            numArrayProteins = ngtmDB.count()
        elif tissue == 'intersection':
            numTissues = len(FUNCTIONAL_TISSUE_LIST)
            numArrayProteins = ngtmDB.find( { 'tissue_list' : 
                                              { '$size' : 
                                                numTissues } } ).count()
        else:
            numArrayProteins = ngtmDB.find( { 'tissue_list' : tissue } ).count()

        # Get Num Proteins According to PPI
        tissueSubgraph = getTissueSubgraph(tissue)
        numPPIProteins = tissueSubgraph.num_vertices()
        
        coverage = float(numPPIProteins)/float(numArrayProteins)

        outFile.write('%s\t%d\t%d\t%f\n' % (tissue, numArrayProteins, 
                                        numPPIProteins, coverage))

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
    
