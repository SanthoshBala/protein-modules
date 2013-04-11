#! /usr/bin/python

# Santhosh Balasubramanian
# March 6, 2013

from settings import *
from pymongo import *

# getEntrezFromEnsembl
def getEntrezFromEnsembl(ensembl_gene_id):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    gim = db.geneIDMap

    record = gim.find_one( { 'ensembl_gene_id' : ensembl_gene_id } )
    if record == None:
        return None
    entrezID = record.get('entrez_gene_id')
    if entrezID == None:
        return None
    else:
        return entrezID[0]


# createModuleAnnotationDB
def createModuleAnnotationDB():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    maDB = db.moduleAnnotation
    gaDB = db.geneAnnotation

    # Iterate through Module DB
    for record in modDB.find():
        moduleID = record.get('module_id')
        print moduleID
        geneList = record.get('gene_list')

        newModuleRecord = {}
        newModuleRecord.update( { 'module_id' : moduleID } )

        moduleFunctions = set()
        moduleProcesses = set()
        moduleComponents = set()

        # Iterate through Genes in this Module
        for geneID in geneList:
            # Get Entrez Gene ID
            if 'ENSG' in geneID:
                geneID = getEntrezFromEnsembl(geneID)
                if geneID == None:
                    continue
            
            # Get Gene Annotations
            annotRecord = gaDB.find_one( { 'entrez_gene_id' : geneID } )
            if annotRecord == None:
                continue
            
            geneFunctionSet = set()
            geneProcessSet = set()
            geneComponentSet = set()

            geneFunctionList = annotRecord.get('function')
            geneProcessList = annotRecord.get('process')
            geneComponentList = annotRecord.get('component')

            for function in geneFunctionList:
                geneFunctionSet.add(tuple(function))
            for process in geneProcessList:
                geneProcessSet.add(tuple(process))
            for component in geneComponentList:
                geneComponentSet.add(tuple(component))

            moduleFunctions = moduleFunctions.union(geneFunctionSet)
            moduleProcesses = moduleProcesses.union(geneProcessSet)
            moduleComponents = moduleComponents.union(geneComponentSet)
            
        newModuleRecord.update( { 'component' : list(moduleComponents) } )
        newModuleRecord.update( { 'function' : list(moduleFunctions) } )
        newModuleRecord.update( { 'process' : list(moduleProcesses ) } )

        maDB.save(newModuleRecord)
    
    cli.close()

# Get Annotation Averages
def getAnnotationAverages():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    maDB = db.moduleAnnotation

    geneSum = 0
    functionSum = 0
    componentSum = 0
    processSum = 0


    for record in maDB.find():
#        geneSum = geneSum + len(record.get('gene_list'))
        functionSum = functionSum + len(record.get('function'))
        processSum = processSum + len(record.get('process'))
        componentSum = componentSum + len(record.get('component'))

    numModules = maDB.count()
#    geneAvg = geneSum/float(numModules)
    functionAvg = functionSum/float(numModules)
    processAvg = processSum/float(numModules)
    componentAvg = componentSum/float(numModules)

#    print geneAvg
    print functionAvg
    print processAvg
    print componentAvg

    cli.close()
