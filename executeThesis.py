#! /usr/bin/python

# executeThesis.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: March 24, 2013

# Python Imports
import sys
import time

# Global Imports
from settings import *

# Utility Imports
from common.shell import *
from common.databases import *
from common.parallelization import *

# Nomenclature Imports
from nomenclature.geneID import *

# Expression Imports
from expression.probeGeneMap import *
from expression.sampleProbeMap import *
from expression.tissueProbeMap import *
from expression.tissueGeneMap import *
from expression.geneTissueMap import *

# Graph Imports

# Networks Imports
from networks.networkCreation import *

# Modules Imports
from modules.moduleUtil import *
from modules.moduleTopology import *

# Analysis Imports
from analysis.proteomeCoverage import *

# analyzeProteomeCoverage: Generates files listing percentage of known
# proteome that is represented in each tissue network.
def analyzeProteomeCoverage(outFile):
    
    # Normal Networks
    startTime = time.time()
    outFile.write('Analyzing Proteome Coverage...')

    getProteomeCoverage()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Shuffled Networks
    startTime = time.time()
    outFile.write('Analyzing Shuffled Proteome Coverage...')

    getProteomeCoverage(shuffle=True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
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

    # Create Global, Intersection, and Tissue Networks
    createGlobalAndTissueNetworks(outFile)

    # Shuffle Gene Expression Data
    shuffleGeneExpression(outFile)

    # Create Shuffled Global, Intersection, and Tissue Networks
    createShuffledNetworks(outFile)
    
    # Analyze Proteome Coverage
    analyzeProteomeCoverage(outFile)

    # Construct Modules for All Tissue Subgraphs
    startTime = time.time()
    outFile.write('Finding Modules for Tissue Subgraphs...\n')
#    findModulesForAllTissues()
    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('Time Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.write('Success: Found Modules\n\n')
    outFile.flush()

    # Construct Module Topology for All Tissue Subgraphs
    startTime = time.time()
    outFile.write('Generating Module Topologies...\n')
#    constructAllModuleTopologies(outFile)
    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('Time Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.write('Success: Generated Tissue Subgraphs\n\n')
    outFile.flush()

    
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


# createAndAnalyzeNomenclatureDBs: Creates and analyzes nomenclature DBs.
def createAndAnalyzeNomenclatureDBs(outFile):

    startTime = time.time()
    outFile.write('Creating Nomenclature DBs...')

    # Drop DBs if They Exist
    dropNomenclatureDBs()

    # Create DB
    createGeneIDMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    startTime = time.time()
    outFile.write('Analyzing ID Map Coverage...')

    # Analyze DB
    analyzeIDMapCoverage()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    return

# createAndAnalyzeProbeGeneMap: Executes necessary helper functions to create
# MongoDB databases that map microarray probes to genes in various 
# nomenclatures, as well as analyzing them to determine the coverage of
# each nomenclature.
def createAndAnalyzeProbeGeneMap(outFile):
    # Drop DBs if Exist
    dropProbeGeneMapDBs()
    
    # Go from Microarray Annotation File to
    # { probe : { ensemble: id, unigene: id, entrez : id,...} }
    startTime = time.time()
    outFile.write('Creating Probe:Gene Map...')
    
    # Create Probe:Gene Map
    createProbeGeneMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    startTime = time.time()
    outFile.write('Analyzing Probe:Gene Map Coverage...')

    # Analyze Probe Gene Map Coverage
    analyzeProbeGeneMapCoverage()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Extend Probe:Gene Mapping
    # Combines Probe:Gene Map with GeneID Map
    startTime = time.time()
    outFile.write('Creating Extended Probe:Gene Map...\n')
    outFile.write('\tWarning: Will take ~2 hours.\n')

    # Find # of Records in <geneIDMap> for Parallelization
    numGeneRecords = MongoClient().db.geneIDMap.count()
    parallelizeTaskByRecord(createExtProbeGeneMap, numGeneRecords)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()    

    # Merge Extra Records of <probeGeneMap> with <extProbeGeneMap>
    startTime = time.time()
    outFile.write('Merging PGM with EPGM...')

    # Find # of Records in <probeGeneMap> for Parallelization
    numProbeRecords = MongoClient().db.probeGeneMap.count()
    parallelizeTaskByRecord(mergePGMwithEPGM, numProbeRecords)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Analyze <extProbeGeneMap>
    startTime = time.time()
    outFile.write('Analyzing Extended Probe:Gene Map Coverage...')

    analyzeProbeGeneMapCoverage(extended=True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))

    outFile.flush()

    return

# createAndAnalyzeTissueProbeMap: Executes all helper functions necessary
# to construct DB that contains the expression E of probe P in tissue T.
# Analyzes this DB to determine # of expressed, unexpressed, and ambiguous
# probes within each tissue.
def createAndAnalyzeTissueProbeMap(outFile):

    dropTissueProbeMapDBs()

    # Go from Microarray SOFT file to
    # { sample : { probe : [ exp1, exp 2],... }, ... }
    startTime = time.time()
    outFile.write('Generating Sample:Probe Map...')

    # Create Sample:Probe Map
    createSampleProbeMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()

    # Go from { sample : { probe : [exp1, exp2],... },... } to
    #         { tissue : { probe : [exp1, exp2],... },... }
    startTime = time.time()
    outFile.write('Creating Tissue:Probe Map...')

    # Create Tissue:Probe Map
    parallelizeTaskByTissue(createTissueProbeMap)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    startTime = time.time()
    outFile.write('Analyzing Tissue:Probe Map...\n')

    # Analyze Tissue:Probe Map
    analyzeTissueProbeMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('Time Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.write('Success: Analyzed Tissue:Probe Map\n\n')
    outFile.flush()

    # Normalize Tissue:Probe Mapping
    # Goes from { tissue : { probe : [exp1, exp2],... },... } to
    #           { tissue : { probe : expVal,... },... }
    startTime = time.time()
    outFile.write('Normalizing Tissue:Probe Map...\n')

    # Normalize Tissue:Probe Map
    normalizeTissueProbeMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('Time Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.write('Success: Normalized Tissue:Probe Map\n\n')
    outFile.flush()

    return

# createAndAnalyzeGeneTissueMap: Executes all helper functions necessary to
# generate <tissueGeneMap>, <normTissueGeneMap>, <geneTissueMap>, and 
# <normGeneTissueMap>. <normGeneTissueMap> is the only DB that deals with
# presence/absence instead of raw expression values.
def createAndAnalyzeGeneTissueMap(outFile):
    
    dropGeneTissueMapDBs()

    # Generate Tissue:Gene Map
    # Go from { tissue : { 123_at : 123, 234_at : 234, ...}, ... } to
    #         { tissue : { 'entrez_gene' : { gene_id : [ 123, 234], }, }, }
    startTime = time.time()
    outFile.write('Creating Tissue:Gene Map...\n')
    outFile.write('\tWarning: Will take ~10 hours. Printing progress...\n')
    
    parallelizeTaskByTissue(createTissueGeneMap)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.write('DONE\n\n')
    outFile.flush()    

    startTime = time.time()
    outFile.write('Analyzing Tissue:Gene Map...')

    analyzeTissueGeneMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')    
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    
    outFile.flush()

    # Normalize Tissue:Gene Map
    # Go from { tissue : { nomenclature : { gene : [1, 2, 3] } } } to
    #         { tissue : { nomenclature : { gene : 2 } } }
    startTime = time.time()
    outFile.write('Normalizing Tissue:Gene Map...')

    parallelizeTaskByTissue(normalizeTissueGeneMap)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))

    outFile.flush()

    # Create Gene:Tissue Map
    # Go from { tissue : { nomenclature : { gene : 2 } } } to
    #         { gene : { tissue : 2 } }
    startTime = time.time()
    outFile.write('Creating Gene:Tissue Map...')

    record = MongoClient().db.normTissueGeneMap.find_one()
    numGenes = len(record.get(GENE_TISSUE_MAP_NOMENCLATURE).keys())
    parallelizeTaskByItem(createGeneTissueMap, numGenes)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')    
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    
    outFile.flush()

    # Normalize Gene:Tissue Map
    # Go from { gene : { tissue : 2 } } to
    #         { gene : id, tissue_list: [ A, B, C ], }
    startTime = time.time()
    outFile.write('Normalizing Gene:Tissue Map...')

    normalizeGeneTissueMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')    
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    
    outFile.flush()


    # Shuffle Gene:Tissue Map
    # Randomizes gene expression
    startTime = time.time()
    outFile.write('Shuffling Gene:Tissue Map...')

    shuffleGeneTissueMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')    
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    
    outFile.flush()

    return

# createGlobalAndTissueNetworks: Executes helper functions to create
# global, intersection, and tissue networks.
def createGlobalAndTissueNetworks(outFile):
    # Create Global Network
    startTime = time.time()
    outFile.write('Creating Global Network...')

    createGlobalNetwork()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Intersection Network
    startTime = time.time()
    outFile.write('Creating Intersection Network...')

    createIntersectionNetwork()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Tissue Networks
    startTime = time.time()
    outFile.write('Creating Tissue Networks...')

    createTissueNetworks()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    return

# shuffleGeneExpressionData: Shuffles expression data in <normGeneTissueMap>
# and works backward to create <shuffleGeneTissueMap> and 
# <shuffleNormTissueGeneMap>.
def shuffleGeneExpressionData(outFile):

    # Shuffle Gene Tissue Map
    startTime = time.time()
    outFile.write('Shuffling Gene Tissue Map...')

    shuffleGeneTissueMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Normalize <shuffleGeneTissueMap>
    startTime = time.time()
    outFile.write('Normalizing Shuffled Gene Tissue Map...')

    normalizeGeneTissueMap(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Shuffled Normalized Tissue Gene Map
    startTime = time.time()
    outFile.write('Shuffling Normalized Tissue Gene Map...')

    shuffleNormTissueGeneMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    return


# createShuffledNetworks: Executes helper functions necessary to
# create tissue network files.
def createShuffledNetworks(outFile):

    # Create Global Networks
    startTime = time.time()
    outFile.write('Creating Shuffled Global Network...')

    createGlobalNetwork(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Intersection Network
    startTime = time.time()
    outFile.write('Creating Shuffled Intersection Network...')

    createIntersectionNetwork(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Tissue Network
    startTime = time.time()
    outFile.write('Creating Shuffled Tissue Networks...')

    createTissueNetworks(shuffleVal = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('Time Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()
    
    return



# - - - - - - - - - - MAIN SCRIPT - - - - - - - - - - #


def main():
    f = open('testfile7', 'w')
    createGlobalAndTissueNetworks(f)
    shuffleGeneExpressionData(f)
    createShuffledNetworks(f)
    analyzeProteomeCoverage(f)
    return

main()

