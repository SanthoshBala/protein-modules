#! /usr/bin/python

# geneOntology.py
# Author: Santhosh Balasubramanian
# Created: March 6, 2013
# Last Modified: April 16, 2013

# Python Imports
import sys

# Library Imports
from fisher import *
from pymongo import *
from Bio import Entrez

# Global Imports
from settings import *

# Utility Imports
from common.strings import *

# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

# *Always* tell NCBI who you are
Entrez.email = 'user@email.com'

# Inferred from Physical Interaction, Inferred from Direct Assay,
# Traceable Author Statement, Inferred by Curator, Inferred from Experiment
ACCEPTABLE_GO_EVIDENCE = ['IPI', 'IDA', 'TAS', 'IC', 'EXP']

QUERY_GO_TYPES = [ 'function', 'process' ] 

GENE_ANNOTATION_FILENAME = 'human.gene2go'
PATH_TO_GENE_ONTOLOGY = PATH_TO_DATA + 'ontology/'

# - - - - - - - - - - API CALL - - - - - - - - - - #

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

# - - - - - - - - - - DATABASE - - - - - - - - - - #

# parseNCBIGOFile: Extract Human gene annotations from NCBI file.
def parseNCBIGOFile():
    inFile = open(PATH_TO_ONTOLOGY + RAW_GENE_ONTOLOGY_FILENAME, 'r')
    outFile = open(PATH_TO_ONTOLOGY + GENE_ONTOLOGY_FILENAME, 'w')

    for line in inFile:
        # Keep Format Line
        if line[0] == '#':
            outFile.write(line)

        # Get Human Lines Only
        lineFields = parseTabSeparatedLine(line)
        if lineFields[0] == '9606':
            outFile.write(line)
        
    # Close Files
    inFile.close()
    outFile.close()

    return
    
# createGeneOntologyDB
def createGeneOntologyDB():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    coll = db.geneOntology

    # Open Input File
    inFile = open(PATH_TO_ANNOTATION + GENE_ONTOLOGY_FILENAME, 'r')

    # Track Entrez ID as Iterating
    currentEntrezID = '1'
    record = { 'entrez_gene_id' : '1',
               'function' : [],
               'component' : [],
               'process' : [],
        }

    # Iterate through <inFile>
    for line in inFile:
        if line[0] == '#':
            continue
    
        lineFields = parseTabSeparatedLine(line)

        entrezID = lineFields[1]
        goID = lineFields[2]
        goEvidence = lineFields[3]
        goType = lineFields[7].lower()
        # Append Modifier to GO Term
        if lineFields[4] == '-':
            goTerm = lineFields[5]
        else:
            goTerm = lineFields[4] + '_' + lineFields[5]

        # Save Record to DB
        if entrezID != currentEntrezID:
           coll.insert(record)
           record = { 'entrez_gene_id' : entrezID,
                      'function' : [],
                      'component' : [],
                      'process' : [],
                      }
           currentEntrezID = entrezID
           
        # Reject Weak Evidence
        if goEvidence not in ACCEPTABLE_GO_EVIDENCE:
            continue
                    
        if (goID, goTerm) not in record[goType]:
            record[goType].append( (goID, goTerm) )
        
    # Save Last Record
    coll.insert(record)

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return

# getGeneAnnotation: Return all annotations in a single list.
def getGeneAnnotation(geneID, nomenclature = 'entrez', combined = True):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneOntology = db.geneOntology

    # Get Record
    if nomenclature == 'entrez':
        record = geneOntology.find_one( { 'entrez_gene_id' : geneID } )
    elif nomenclature == 'ensembl':
        entrezID = getEntrezForEnsemblGene(geneID)
        if not entrezID:
            return list()
        record = geneOntology.find_one( { 'entrez_gene_id' : entrezID } )
        
    if not record:
        return list()

    functionList = record.get('function')
    processList = record.get('process')
    componentList = record.get('component')

    if combined:
        returnList = []
        if 'function' in QUERY_GO_TYPES:
            returnList = returnList + functionList
        if 'process' in QUERY_GO_TYPES:
            returnList = returnList + processList
        if 'component' in QUERY_GO_TYPES:
            returnList = returnList + componentList
        return returnList
    else:
        numTypes = len(QUERY_GO_TYPES)
        if numTypes == 1:
            return functionList
        if numTypes == 2:
            return functionList, processList
        if numTypes == 3:
            return functionList, processList, componentList

