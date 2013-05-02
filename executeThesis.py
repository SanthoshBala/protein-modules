#! /usr/bin/python

# executeThesis.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: March 24, 2013

# Global Imports
from settings import *

# Helper Functions
from thesis import *


# createAndAnalyzeDeBruijnGraph
def createAndAnalyzeDeBruijnGraph(outFile):


    # Create De Bruijn Graph
    startTime = time.time()
    outFile.write('Creating De Bruijn Graph...\n')

    tissueList = FUNCTIONAL_TISSUE_LIST + ['intersection']
    createDeBruijnSetGraph(tissueList = tissueList)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()

    # Create Edge Matrices
    startTime = time.time()
    outFile.write('Creating Edge Matrices...\n')

    getEmbryonicDeBruijnAdditionMatrix()
    getEmbryonicDeBruijnRemovalMatrix()
    getEmbryonicDeBruijnExchangeMatrix()
    getDeBruijnEdgeTypeCounts()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()


    return


# executeThesis: Executes all helper functions necessary to create databases 
# and data files for tissue-specific analysis of PPI networks.
def executeThesis():

    # Open Output File
    outFile = open(PATH_TO_ANALYSIS + OUTPUT_FILENAME, 'w')

    # Add Working Directory to PYTHONPATH
    outFile.write('Mounting code directory...')

    mountCodeDir()

    outFile.write('DONE\n\n')
    outFile.flush()

    # Create/Analyze Gene Nomenclature DBs
    createAndAnalyzeNomenclatureDBs(outFile)

    # Create/Analyze Probe:Gene Map
    createAndAnalyzeProbeGeneMap(outFile)

    # Create/Analyze Tissue:Probe Map
    createAndAnalyzeTissueProbeMap(outFile)

    # Create/Analyze Gene:Tissue Map
    createAndAnalyzeGeneTissueMap(outFile)

    # Shuffle Gene Expression Data
    shuffleGeneExpression(outFile)

    # Create Global, Intersection, and Tissue Networks
    createGlobalAndTissueNetworks(outFile)

    # Create Shuffled Global, Intersection, and Tissue Networks
    createShuffledNetworks(outFile)
    
    # Find SPICI Modules For All Networks
    findRawModulesForAllNetworks(outFile)

    # Create Module Database
    createModuleDatabase(outFile)
    createModuleIDFilesForAllNetworks(outFile)

    # Analyze Proteome Coverage
    analyzeProteomeCoverage(outFile)

    # Analyze Tissue Multiplicity
    analyzeTissueMultiplicity(f)

    # Create and Analyze Module Topologies
    createAndAnalyzeModuleTopologies(outFile)

    # Create and Analyze De Bruijn Graph


    
    # Construct Housekeeping Genes DB
    constructHousekeepingGenesDB()
    fillHousekeepingGenesDB()
    fillHousekeepingGeneModules()
    writeHousekeepingGeneModuleIdentitiesToFile()
    writeHousekeepingGeneModuleGermLayersToFile()


    # Close output and error files
    outFile.close()
    errFile.close()


# - - - - - - - - - - HELPER FUNCTIONS - - - - - - - - - - #















# - - - - - - - - - - MAIN SCRIPT - - - - - - - - - - #


def main():
    f = open('testfile9', 'w')
    createAndAnalyzeDeBruijnGraph(f)
    return

main()
