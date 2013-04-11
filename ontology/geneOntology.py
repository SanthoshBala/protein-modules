#! /usr/bin/python

# Santhosh Balasubramanian
# March 6, 2013

import sys

from settings import *
from common import *
from pymongo import *
from Bio import Entrez

ACCEPTABLE_GO_EVIDENCE = ['IPI', 'IDA', 'TAS', 'IC']

# *Always* tell NCBI who you are
Entrez.email = 'santhosh@princeton.edu'

# Call NCBI Entrez DB to get gene annotation for gene (by Entrez ID).
def getNCBIGeneAnnotation(entrez_gene_id):
    """ Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information. Returns
    a list of dictionaries with the annotations."""
    
    idList = [entrez_gene_id]

    request = Entrez.epost("gene", id= ','.join(idList))
    
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print 'An error occurred while retrieving the annotations.'
        print 'The error returned was %s' % e
        sys.exit(-1)

    webEnv = result['WebEnv']
    queryKey = result['QueryKey']
    data = Entrez.esummary(db='gene', webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    return annotations

# createAnnotationDB
def createAnnotationDB():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    coll = db.geneAnnotation

    # Open Input File
    inFile = open(PATH_TO_ANNOTATION + 'human.gene2go', 'r')
    
    currentEntrezID = '1'
    record = { 'entrez_gene_id' : '1',
               'function' : [],
               'component' : [],
               'process' : [],
        }

    for line in inFile:
        # Ignore Comments
        if line[0] == '#':
            continue
    
        lineFields = parseTabSeparatedLine(line)

        entrezID = lineFields[1]

        if entrezID != currentEntrezID:
           coll.insert(record)
           record = { 'entrez_gene_id' : entrezID,
                      'function' : [],
                      'component' : [],
                      'process' : [],
                      }
           currentEntrezID = entrezID
           
        goID = lineFields[2]
        
        goEvidence = lineFields[3]
        if goEvidence not in ACCEPTABLE_GO_EVIDENCE:
            continue

        if lineFields[4] == '-':
            goTerm = lineFields[5]
        else:
            goTerm = lineFields[4] + '_' + lineFields[5]
            
        goType = lineFields[7].lower()

        if (goID, goTerm) not in record[goType]:
            record[goType].append( (goID, goTerm) )
        

# getEntrezGeneAnnotation
def getEntrezGeneAnnotation(entrez_gene_id):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneAnnotation = db.geneAnnotation

    # Get Record
    record = geneAnnotation.find_one( { 'entrez_gene_id' : entrez_gene_id } )
    if not record:
        return list()

    functionList = record.get('function')
    processList = record.get('process')

    return functionList + processList

# getEntrezGeneFunctionsProcesses
def getEntrezGeneFunctionsProcesses(entrez_gene_id):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneAnnotation = db.geneAnnotation

    # Get Record
    record = geneAnnotation.find_one( { 'entrez_gene_id' : entrez_gene_id } )
    if not record:
        return None, None

    functionList = record.get('function')
    processList = record.get('process')

    return functionList, processList


# getEnsemblGeneAnnotation
def getEnsemblGeneAnnotation(ensembl_gene_id):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneAnnotation = db.geneAnnotation
    geneIDMap = db.geneIDMap

    # Get Entrez ID
    geneRecord = geneIDMap.find_one( { 'ensembl_gene_id' : ensembl_gene_id } )

    if not geneRecord:
        return None, None

    entrezID = geneRecord.get('entrez_gene_id')
    if not entrezID:
        return None, None
    entrezID = entrezID[0]

    # Get Record
    record = geneAnnotation.find_one( { 'entrez_gene_id' : entrezID } )
    if not record:
        return None, None

    functionList = record.get('function')
    processList = record.get('process')

    return functionList + processList

# getEnsemblGeneFunctionsProcesses
def getEnsemblGeneFunctionsProcesses(ensembl_gene_id):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneAnnotation = db.geneAnnotation
    geneIDMap = db.geneIDMap

    # Get Entrez ID
    geneRecord = geneIDMap.find_one( { 'ensembl_gene_id' : ensembl_gene_id } )

    if geneRecord == None:
        return None, None

    entrezID = geneRecord.get('entrez_gene_id')
    if entrezID == None:
        return None, None
    entrezID = entrezID[0]

    # Get Record
    record = geneAnnotation.find_one( { 'entrez_gene_id' : entrezID } )
    if not record:
        return None, None

    functionList = record.get('function')
    processList = record.get('process')

    return functionList, processList

    
