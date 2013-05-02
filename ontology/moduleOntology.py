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

# Ontology Imports
from ontology.geneOntology import *

# createModuleOntologyDB
def createModuleOntologyDB(skipVal = 0, limitVal = 1000):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules
    modOntDB = db.moduleOntology
    geneOntDB = db.geneOntology

    # Iterate through Module DB
    for record in modDB.find(skip=skipVal, limit=limitVal):
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
        
        numModGenes = len(entrezGeneList)
        numGenes = geneOntDB.find().count()

        # Iterate through Genes in <moduleID>
        for geneID in entrezGeneList:
            # Get Gene Annotations
            geneFuncs, geneProcs, geneComps = getGeneAnnotation(geneID,
                                                                combined = False)
            
            # Iterate through Function Annotations
            for funcAnnot in geneFuncs:
                # Fill in Matrix for Fisher's Exact Test                
                query = constructQueryForGO(funcAnnot, 'function',
                                            entrezGeneList)
                # Number of Genes in Module with Annotation
                numModPos = geneOntDB.find(query).count()
                # Number of Genes in Module without Annotation
                numModNeg = numModGenes - numModPos

                query = constructQueryForGO(funcAnnot, 'function')
                # Number of Genes not in Module with Annotation
                numGenPos = geneOntDB.find(query).count() - numModPos
                # Number of Genes not in Module without Annotation
                numGenNeg = numGenes - numGenPos - numModPos

                p = pvalue(numModPos, numModNeg, numGenPos, numGenNeg)

                if p.two_tail < 0.05:
                    modFunctions.add(tuple(funcAnnot))
    
            # Iterate through Process Annotations
            for procAnnot in geneProcs:
                # Fill in Matrix for Fisher's Exact Test                
                query = constructQueryForGO(procAnnot, 'process',
                                            entrezGeneList)
                # Number of Genes in Module with Annotation
                numModPos = geneOntDB.find(query).count()
                # Number of Genes in Module without Annotation
                numModNeg = numModGenes - numModPos

                query = constructQueryForGO(procAnnot, 'process')
                # Number of Genes not in Module with Annotation
                numGenPos = geneOntDB.find(query).count() - numModPos
                # Number of Genes not in Module without Annotation
                numGenNeg = numGenes - numGenPos - numModPos

                p = pvalue(numModPos, numModNeg, numGenPos, numGenNeg)

                if p.two_tail < 0.05:
                    modProcesses.add(tuple(procAnnot))

            # Iterate through Component Annotations
            for compAnnot in geneComps:
                # Fill in Matrix for Fisher's Exact Test                
                query = constructQueryForGO(compAnnot, 'component', 
                                            entrezGeneList)
                # Number of Genes in Module with Annotation
                numModPos = geneOntDB.find(query).count()
                # Number of Genes in Module without Annotation
                numModNeg = numModGenes - numModPos

                query = constructQueryForGO(compAnnot, 'component')
                # Number of Genes not in Module with Annotation
                numGenPos = geneOntDB.find(query).count() - numModPos
                # Number of Genes not in Module without Annotation
                numGenNeg = numGenes - numGenPos - numModPos

                p = pvalue(numModPos, numModNeg, numGenPos, numGenNeg)

                if p.two_tail < 0.05:
                    modComponents.add(tuple(compAnnot)) 
                    
        modOntologyRecord.update( { 'component' : list(modComponents) } )
        modOntologyRecord.update( { 'function' : list(modFunctions) } ) 
        modOntologyRecord.update( { 'process' : list(modProcesses) } )
        modOntDB.save(modOntologyRecord)
    
    # Close DB Connection    
    cli.close()
    
    return

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
