#! /usr/bin/python

# moduleOntology.py
# Author: Santhosh Balasubramanian
# Created: March 6, 2013
# Last Modified: April 16, 2013


# Library Imports
from fisher import *
from pymongo import *

# Global Imports
from settings import *

# Nomenclature Imports
from nomenclature.geneID import *


# createModuleOntologyDB
def createModuleOntologyDB(hypergeometric = True):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    modOntDB = db.moduleOntology
    geneOntDB = db.geneOntology

    # Iterate through Module DB
    for record in modDB.find():
        moduleID = record.get('module_id')
        geneList = record.get('gene_list')
        if PRINT_PROGRESS:
            print moduleID
        
        # Initialize <moduleOntology> Record
        modOntologyRecord = { 'module_id' : moduleID }
        modFunctions = set()
        modProcesses = set()
        modComponents = set()
        rejectedAnnotations = set()

        # Get Entrez Gene List
        entrezGeneList = []
        for geneID in geneList:
            if 'ENSG' in geneID:
                entrezID = getEntrezForEnsemblGene(geneID)
                if not entrezID:
                    continue
            else:
                entrezID = geneID
            entrezGeneList.append(entrezID)
        
        # Iterate through Genes in <moduleID>
        for geneID in entrezGeneList:
            # Get Gene Annotations
            geneFuncs, geneProcs, geneComps = getGeneAnnotation(geneID,
                                                                combined = False)
            
            for funcAnnot in geneFuncs:
                # Fill in Matrix for Fisher's Exact Test                
                constructQueryForGO('function', funcAnnot, entrezGeneList)

        # Iterate through Genes in <moduleID>
        for geneID in geneList:
            # Get Gene Annotations
            geneOntologyRecord = goDB.find_one( { 'entrez_gene_id' : geneID } )
            if not geneOntologyRecord:
                continue

            geneFunctions = set()
            geneProcesses = set()
            geneComponents = set()

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
