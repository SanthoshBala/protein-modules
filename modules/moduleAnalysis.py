#! /usr/bin/python/

# moduleAnalysis.py
# Author: Santhosh Balasubramanian
# Date: February 1, 2013

from settings import *
from common import *
from pymongo import *
from random import shuffle
import subprocess
import os
from graphs.graphIO import *
from ontology.geneOntology import *
from graphs.graphUtil import *
from modules.moduleUtil import *

JACCARD_THRESHOLD = 0.85
JACCARD_INT = 85

def getAllModuleSimilarities():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    mtm = db.moduleTissueMap

    # Open Outfile
    outFilename = PATH_TO_ANALYSIS + 'similarity.all.module.list'
    outFile = open(outFilename, 'w')

    for recordA in mtm.find():
        moduleA = recordA.get('module_id')
        print moduleA
        setA = set(recordA.get('gene_list'))
        
        for recordB in mtm.find():
            moduleB = recordB.get('module_id')
            if moduleB < moduleA:
                continue
            if moduleA == moduleB:
                continue
            
            setB = set(recordB.get('gene_list'))
            
            jaccard = getJaccardSimilarity(setA, setB)
            
            outFile.write('%f\n' % jaccard)

    outFile.close()
    cli.close()



def findSimilarModules(jaccardDecimal, jaccardInt):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    mtm = db.moduleTissueMap

    # Open Outfile
    
    outFilename = PATH_TO_ANALYSIS + 'similarity.%s.module.list' % (jaccardInt)
    outFile = open(outFilename, 'w')

    for recordA in mtm.find():
        moduleA = recordA.get('module_id')
        print moduleA
        setA = set(recordA.get('gene_list'))
        
        for recordB in mtm.find():
            moduleB = recordB.get('module_id')
#            print '\t%s' % moduleB
            if moduleB < moduleA:
                continue
            if moduleA == moduleB:
                continue
            
            setB = set(recordB.get('gene_list'))
            
            jaccard = getJaccardSimilarity(setA, setB)
            
            if jaccard > jaccardDecimal:
                outFile.write('%s\t%s\n' % (moduleA, moduleB))

    outFile.close()
    cli.close()


# Get Module-Tissue Histogram
def getModuleTissueHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    mtm = db.modules

    # Open Outfile
    outFile = open(PATH_TO_ANALYSIS + 'module.tissue.histogram', 'w')

    histList = [0]*(len(CANONICAL_TISSUE_LIST) + 1)

    # Iterate through module tissue map
    for record in mtm.find():
        numTissues = len(record.get('tissue_list'))
        histList[numTissues] = histList[numTissues] + 1

    for i in range(len(histList)):
        outFile.write('%d\t%d\n' % (i, histList[i]))

# getModuleOrganSystemHistogram
def getModuleOrganSystemHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Outfile
    outFile = open(PATH_TO_ANALYSIS + 'module.organ_system.histogram', 'w')
    
    histList = [0]*(len(CANONICAL_ORGAN_SYSTEM_LIST) + 1)

    # Iterate through modules DB
    for record in modDB.find():
        organSet = set()
        tissueList = record.get('tissue_list')
        for tissue in tissueList:
            organ = TISSUE_SYSTEM_MAP.get(tissue)
            organSet.add(organ)
        numOrgans = len(organSet)
        histList[numOrgans] = histList[numOrgans] + 1

    for i in range(len(histList)):
        outFile.write('%d\t%d\n' % (i, histList[i]))

    outFile.close()


# analyzeAllSimilarityTopologies
def analyzeAllSimilarityTopologies():
    for i in range(5, 100, 5):
        print i
        if i == 5:
            analyzeModuleSimilarityTopology('05')
        else:
            analyzeModuleSimilarityTopology(str(i))

# analyzeModuleSimilarityTopology
def analyzeModuleSimilarityTopology(threshold):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    maDB = db.moduleAnnotation

    # Open Input File
    inFilename = MODULE_SIMILARITY_BASE_FILENAME % threshold
    inFilepath = PATH_TO_MODULE_SIMILARITY + inFilename
        
    # Open Output File
    outFilename = 'similarity.%s.topology.analysis' % threshold
    outFilepath = PATH_TO_MODULE_SIMILARITY + outFilename
    outFile = open(outFilepath, 'w')

    # Read Similarity Graph
    simGraph = readCanonicalGraph(inFilepath)

    # Get Connected Components
    componentProp, histogram = getConnectedComponents(simGraph)
    numConnectedComponents = len(histogram)

    # Iterate Through each Connected Component
    geneIDProp = simGraph.vertex_properties['gene_id']
    for i in range(numConnectedComponents):
        componentVertices = getVertexSetByProperty(simGraph, componentProp, i)
        
        outFile.write('COMPONENT %d\n' % i)
        
        # Iterate through Vertices of Component
        for vertex in componentVertices:
            moduleID = geneIDProp[vertex]
            
            outFile.write('\t%s\n' % moduleID)
            
            # Get Annotations
            record = maDB.find_one( { 'module_id' : moduleID } )
            functions = record.get('function')
            processes = record.get('process')
            annotations = functions + processes
            
            # Print Annotations
            for annotation in annotations:
                outFile.write('\t\t%s\n' % str(annotation))

        outFile.write('\n')

# - - - - - Module Identity - - - - - #

# updateHousekeepingGenesModules
def fillHousekeepingGeneModules():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    hgDB = db.housekeepingGenes
    modDB = db.modules

    # Iterate through all Housekeeping Genes
    for geneRecord in hgDB.find():
        geneID = geneRecord.get('gene_id')
        
        # Find gene in modDB
        modList = modDB.find( { 'gene_list' : geneID } )
        geneRecord['module_identity'] = {}

        # Add entries in hgDB
        for modRecord in modList:
            moduleID = modRecord.get('module_id')
            tissueList = modRecord.get('tissue_list')
            for tissue in tissueList:
                geneRecord['module_identity'].update( { tissue : moduleID } )

        hgDB.save(geneRecord)

    cli.close()

    return
                
def writeHousekeepingGeneModuleIdentitiesToFile():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    hgDB = db.housekeepingGenes

    # Open Output File
    outFilename = 'housekeeping.gene.module.identities'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    # Iterate through DB
    for record in hgDB.find():
        geneID = record.get('gene_id')
        outFile.write('%s\n' % geneID)
        modIdentities = record.get('module_identity')
        if type(modIdentities) == list:
            continue
        for tissue, modID in modIdentities.iteritems():
            outFile.write('\t%s\t%s\n' % (tissue, modID))
            
    outFile.close()
    return

def writeHousekeepingGeneModuleGermLayersToFile():
# Open DB Connection
    cli = MongoClient()
    db = cli.db
    hgDB = db.housekeepingGenes
    modDB = db.modules

    # Open Output File
    outFilename = 'housekeeping.gene.module.germ_layers'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    # Iterate through DB
    for record in hgDB.find():
        geneID = record.get('gene_id')
        outFile.write('%s\n' % geneID)
        modIdentities = record.get('module_identity')
        if type(modIdentities) == list:
            continue
        for tissue, modID in modIdentities.iteritems():
            modRecord = modDB.find_one( { 'module_id' : modID } )
            germLayers = modRecord.get('germ_layer')
            outFile.write('\t%s\t%s\n' % (modID, str(germLayers)))
            
    outFile.close()
    return
    

# - - - - - GENE MODULE MAP - - - - - #


def writeGeneModuleMapToFile():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    gmmDB = db.geneModuleMap

    # Open Output File
    outFilename = 'gene.module.map.analysis'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    # Iterate through GMM
    for geneRecord in gmmDB.find():
        moduleSet = set()
        for tissue, module in geneRecord.iteritems():
            if tissue == '_id':
                continue
            if tissue == 'gene_id':
                continue

            moduleSet.add(module)

        if len(moduleSet) > 1:
            outFile.write('%s\n' % geneRecord.get('gene_id'))
            for tissue, module in geneRecord.iteritems():
                if tissue == '_id':
                    continue
                if tissue == 'gene_id':
                    continue
                outFile.write('\t%s\t%s\n' % (tissue, module))
                

    outFile.close()
    cli.close()
    return

def getModuleTissueHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilename = 'module.tissue.histogram'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    histogram = {}

    # Iterate through modDB
    for record in modDB.find():
        tissueList = record.get('tissue_list')
        numTissues = len(tissueList)
        if histogram.get(numTissues):
            histogram[numTissues] = histogram[numTissues] + 1
        else:
            histogram.update( { numTissues: 1 } )

    # Write to output
    for bucket, count in histogram.iteritems():
        outFile.write('%d\t%d\n' % (bucket, count))

    outFile.close()
    cli.close()
    return

def getModuleTissueExpression():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilename = 'module.tissue.expression'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    tissueModDict = {}
    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    for tissue in AUGMENTED_TISSUE_LIST:
        tissueModDict.update( { tissue : set() } )

    # Iterate through modDB
    for record in modDB.find():
        moduleID = record.get('module_id')
        tissueList = record.get('tissue_list')
        
        for tissue in tissueList:
            tissueModDict[tissue].add(moduleID)

    # Write to output
    for tissue in AUGMENTED_TISSUE_LIST:
        outFile.write('%s\t%d\n' % (tissue, len(tissueModDict[tissue])))

    outFile.close()
    cli.close()
    return


def getModuleOrganSystemExpression():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilename = 'module.organ_system.expression'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    tissueModDict = {}
    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    for system in AUGMENTED_ORGAN_SYSTEM_LIST:
        tissueModDict.update( { system : set() } )

    # Iterate through modDB
    for record in modDB.find():
        moduleID = record.get('module_id')
        tissueList = record.get('tissue_list')
        
        for tissue in tissueList:
            if tissue == 'global':
                system = 'global'
            else:
                system = TISSUE_SYSTEM_MAP.get(tissue)

            tissueModDict[system].add(moduleID)

    # Write to output
    for tissue in AUGMENTED_ORGAN_SYSTEM_LIST:
        outFile.write('%s\t%d\n' % (tissue, len(tissueModDict[tissue])))

    outFile.close()
    cli.close()
    return

    
def getModuleSizeHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFilename = 'module.size.histogram'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    histogram = {}

    # Iterate through modDB
    for record in modDB.find():
        geneList = record.get('gene_list')
        numGenes = len(geneList)
        if histogram.get(numGenes):
            histogram[numGenes] = histogram[numGenes] + 1
        else:
            histogram.update( { numGenes: 1 } )

    # Write to output
    for bucket, count in histogram.iteritems():
        outFile.write('%d\t%d\n' % (bucket, count))

    outFile.close()
    cli.close()
    return

# writeModuleGermLayerMapToFile
def writeModuleGermLayerMapToFile():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Output File
    outFilename = 'module.germ_layer.attributes'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    for moduleRecord in modDB.find():
        ectoBool = False
        mesoBool = False
        endoBool = False

        moduleID = moduleRecord.get('module_id')
        germLayers = moduleRecord.get('germ_layer')

        if 'ectoderm' in germLayers:
            ectoBool = True
        if 'endoderm' in germLayers:
            endoBool = True
        if 'mesoderm' in germLayers:
            mesoBool = True

        outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (moduleID, 
                                                    endoBool, endoBool, 
                                                    mesoBool, mesoBool, 
                                                    ectoBool, ectoBool))
    
    cli.close()
    outFile.close()



# writeModuleGermLayerMapToFile
def writeModuleGermLayerHashToFile():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Output File
    outFilename = 'module.germ_layer.attributes'
    outFilepath = '/home/santhosh/Dropbox/Thesis/data/analysis/debruijn/' + outFilename
    outFile = open(outFilepath, 'w')

    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        value = getModuleGermLayerHash(moduleID)

        outFile.write('%s\t%d\n' % (moduleID, value)) 

    cli.close()
    outFile.close()


# getGermLayerHashHistogram
def getGermLayerHashHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    
    # Open Output File
    outFilename = 'germ_layer.hash.gt4.histogram'
    outFilepath = PATH_TO_MODULE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    histogram = [0]*8

    for moduleRecord in modDB.find():
        moduleID = moduleRecord.get('module_id')
        value = getGermLayerHash(moduleID)
        
        if len(moduleRecord.get('gene_list')) < 4:
            continue

        histogram[value] = histogram[value] + 1


    # Write to File
    for i in range(8):
        outFile.write('%d\t%d\n' % (i, histogram[i]))
        
    cli.close()
    outFile.close()


