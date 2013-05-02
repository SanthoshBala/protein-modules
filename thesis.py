#! /usr/bin/python

# thesis.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: April 25, 2013


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
from modules.deBruijn import *
from modules.moduleUtil import *
from modules.moduleTopology import *
from modules.moduleDatabase import *

# Analysis Imports
from analysis.proteomeCoverage import *
from analysis.tissueMultiplicity import *
from analysis.topologyAnalysis import *

# Ontology Imports
from ontology.geneOntology import *
from ontology.moduleOntology import *



##############################################################################
##                           NOMENCLATURE DATABASES                         ##
##############################################################################


# createAndAnalyzeNomenclatureDBs: Creates and analyzes nomenclature DBs.
def createAndAnalyzeNomenclatureDBs(outFile):

    startTime = time.time()
    outFile.write('Creating Nomenclature DBs...')

    # Drop DBs if They Exist
    dropNomenclatureDBs()

    # Create DB
    createGeneIDMap()
    createUniprotIDMap()

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
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
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
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
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

    return


# shuffleGeneExpressionData: Shuffles expression data in <normGeneTissueMap>
# and works backward to create <shuffleGeneTissueMap> and 
# <shuffleNormTissueGeneMap>.
def shuffleGeneExpressionData(outFile):

    # Drop Shuffled Expression Databases
    dropShuffleExpressionDBs()

    # Shuffle Gene Tissue Map
    startTime = time.time()
    outFile.write('Shuffling Gene Tissue Map...')

    shuffleGeneTissueMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Normalize <shuffleGeneTissueMap>
    startTime = time.time()
    outFile.write('Normalizing Shuffled Gene Tissue Map...')

    normalizeGeneTissueMap(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Shuffled Normalized Tissue Gene Map
    startTime = time.time()
    outFile.write('Shuffling Normalized Tissue Gene Map...')

    shuffleNormTissueGeneMap()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    return


##############################################################################
##                           NETWORK CREATION                               ##
##############################################################################


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
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Intersection Network
    startTime = time.time()
    outFile.write('Creating Intersection Network...')

    createIntersectionNetwork()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Tissue Networks
    startTime = time.time()
    outFile.write('Creating Tissue Networks...')

    createTissueNetworks()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
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
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Intersection Network
    startTime = time.time()
    outFile.write('Creating Shuffled Intersection Network...')

    createIntersectionNetwork(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Create Tissue Network
    startTime = time.time()
    outFile.write('Creating Shuffled Tissue Networks...')

    createTissueNetworks(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()
    
    return


##############################################################################
##                             PROTEOME COVERAGE                            ##
##############################################################################


# analyzeProteomeCoverage: Generates files listing percentage of known
# proteome that is represented in each tissue network.
def analyzeProteomeCoverage(outFile):
    
    # Normal Networks
    startTime = time.time()
    outFile.write('Analyzing Proteome Coverage...')

    getProteomeCoverage()
    getProteomeCoverage(shuffle=True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Compare Global and Local Networks
    startTime = time.time()
    outFile.write('Comparing Global and Local Networks...')
    
    compareGlobalAndLocalNetworks()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    return


##############################################################################
##                        PROTEIN INTERACTION MODULES                       ##
##############################################################################


# findRawModulesForAllNetworks
def findRawModulesForAllNetworks(outFile):
    # Compute Raw SPICI Module Files
    startTime = time.time()
    outFile.write('Computing Raw SPICI Module Files...')
    
    parallelizeTaskByTissue(createAllRawModuleFiles, augmented = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    # Compute Raw SPICI Shuffled Module Files
    startTime = time.time()
    outFile.write('Computing Shuffle Raw SPICI Module Files...')
    
    parallelizeTaskByTissue(createAllRawModuleFiles, augmented = True, 
                            shuffleParam = True, shuffleVal = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    return


# createModuleDatabase
def createModuleDatabase(outFile):
    
    # Drop Module Database
    dropModuleDBs()
    
    # Compute Raw SPICI Shuffled Module Files
    startTime = time.time()
    outFile.write('Initializing Module DB...')
    
    createModuleDB()
    createModuleDB(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    # Add Germ Layer Field to DB
    startTime = time.time()
    outFile.write('Adding Germ Layer Field to Module DB...')
    
    createModuleGermLayerField()
    createModuleGermLayerField(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Add Protein Universality Field to Module DB
    startTime = time.time()
    outFile.write('Adding Protein Universality Field to Module DB...')
    
    numRecords = MongoClient().db.modules.count()
    parallelizeTaskByRecord(createModuleProteinUniversalityField, numRecords)
    parallelizeTaskByRecord(createModuleProteinUniversalityField, numRecords,
                            shuffleParam = True, shuffleVal = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Add Housekeeping Index Field to Module DB
    startTime = time.time()
    outFile.write('Adding Housekeeping Index Field to Module DB...')
    
    numRecords = MongoClient().db.modules.count()
    parallelizeTaskByRecord(createModuleHousekeepingIndexField, numRecords)
    parallelizeTaskByRecord(createModuleHousekeepingIndexField, numRecords,
                            shuffleParam = True, shuffleVal = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    return

# createModuleIDFIlesForAllNetworks
def createModuleIDFilesForAllNetworks(outFile):

    # Create All Module ID Files
    startTime = time.time()
    outFile.write('Creating All Module ID Files...')
    
    createAllModuleIDFiles()
    createAllModuleIDFiles(shuffle = True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    return


##############################################################################
##                             TISSUE MULTIPLICITY                          ##
##############################################################################

# analyzeTissueMultiplicity: Computes tissue multiplicity for vertices,
# edges, and modules.
def analyzeTissueMultiplicity(outFile):

    # Compute Protein Tissue Multiplicity
    startTime = time.time()
    outFile.write('Computing Protein Tissue Multiplicity...')
    
    getProteinTissueMultiplicity()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Compute Interaction Tissue Multiplicity
    startTime = time.time()
    outFile.write('Computing Interaction Tissue Multiplicity...')
    
    getInteractionTissueMultiplicity()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Compute Module Tissue Multiplicity
    startTime = time.time()
    outFile.write('Computing Module Tissue Multiplicity...')
    
    getModuleTissueMultiplicity()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Compute Module Protein Universality
    startTime = time.time()
    outFile.write('Computing Module Protein Universality...')
    
    getModuleProteinUniversality()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    # Compute Module Housekeeping Index
    startTime = time.time()
    outFile.write('Computing Module Housekeeping Index...')
    
    getModuleHousekeepingIndex()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()

    return


##############################################################################
##                           NETWORK TOPOLOGY                               ##
##############################################################################


# analyzeNetworkTopology: Analyze topologies of all tissue networks.
def analyzeNetworkTopology(outFile):

    # Analyze Mean Path Lengths
    startTime = time.time()
    outFile.write('Analyzing Mean Path Lengths...')
    
    parallelizeTaskByTissue(getNetworkMeanPathLengths, augmented=True)
    
    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    # Analyze Degree Distributions
    startTime = time.time()
    outFile.write('Analyzing Degree Distributions...')
    
    parallelizeTaskByTissue(getNetworkDegreeDistributions, augmented=True)
    
    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    # Analyze LCC vs Degree
    startTime = time.time()
    outFile.write('Analyzing LCC vs Degree...')
    
    parallelizeTaskByTissue(getNetworkDegreeDistributions, augmented=True)
    
    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('DONE\n')
    outFile.write('\tTime Elapsed: %d:%d:%d\n\n' % (h,m,s))
    outFile.flush()
    
    return


# createAndAnalyzeModuleTopologies: Create module topologies, annotate them
# with gene ontologies, and check number of non-leaf modules with annotations.
def createAndAnalyzeModuleTopologies(outFile):

    # Construct Module Topology for All Tissue Subgraphs
    startTime = time.time()
    outFile.write('Generating Module Topologies...\n')

    parallelizeTaskByTissue(createAllModuleTopologies, augmented=True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()

    # Construct Module Topology for All Tissue Subgraphs
    startTime = time.time()
    outFile.write('Analyzing Module Topologies...\n')

    parallelizeTaskByTissue(getModuleTopologyAnnotation, augmented=True)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()
        
    return


# createOntologyDatabases: Creates <geneOntology> and <moduleOntology>.
def createOntologyDatabases(outFile):

    # Drop Ontology DBs
    dropOntologyDBs()

    # Create Gene Ontology DB
    startTime = time.time()
    outFile.write('Creating Gene Ontology DB...\n')

    parseNCBIGOFile()
    createGeneOntologyDB()

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()

    # Create Module Ontology DB
    startTime = time.time()
    outFile.write('Creating Module Ontology DB...\n')

    numRecords = MongoClient().db.modules.count()
    parallelizeTaskByRecord(createModuleOntologyDB, numRecords)

    elapsedTime = time.time() - startTime
    h, m, s = hoursMinutesSeconds(elapsedTime)
    outFile.write('\tTime Elapsed: %d:%d:%d\n' % (h,m,s))
    outFile.flush()

    return
