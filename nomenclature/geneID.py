#! /usr/bin/python/

# geneID.py
# Author: Santhosh Balasubramanian
# Created: January 14, 2013
# Last Modified: March 27, 2013


# Python Imports
import json
import math
from urllib2 import urlopen
from multiprocessing import *

# Library Imports
from pymongo import *
from bs4 import BeautifulSoup

# Global Imports
from settings import *

# Utility Imports
from common.strings import *
from common.databases import *
from common.parallelization import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


## Universal Settings

# Determines which data sources are integrated into <geneIDMap>.


# Maps field names from each datasource to canonical values.
CANONICAL_GENE_FIELDS = {
    # GeneMania Terms
    'Ensembl Gene ID' : 'ensembl_gene_id',
    'Ensembl Protein ID' : 'ensembl_protein_id',
    'Entrez Gene ID' : 'entrez_gene_id',
    'Gene Name' : 'gene_symbol',
    'RefSeq Protein ID' : 'refseq_protein_id',
    'RefSeq mRNA ID' : 'refseq_transcript_id',
    'Synonym' : 'gene_symbol',
    'Uniprot ID' : 'uniprot_id',
    'Primary Gene ID' : 'primary_gene_id',

    # BioGRID Terms
    'ENTREZ_GENE' : 'entrez_gene_id',
    'SYNONYM' : 'gene_symbol',
    'REFSEQ_PROTEIN_ACCESSION' : 'refseq_protein_id',
    'ENSEMBL' : 'ensembl_gene_id', 
    'REFSEQ_DNA_ACCESSION' : 'refseq_gene_id', 
    'OFFICIAL_SYMBOL' : 'gene_symbol', 
    'ENTREZ_GENE' : 'entrez_gene_id', 
    'HGNC' : 'hgnc_id',

    # UniProt Terms
    'RefSeq' : 'refseq_protein_id',
    'Ensembl_PRO' : 'ensembl_protein_id', 
    'Ensembl' : 'ensembl_gene_id',
    'UniGene' : 'unigene_id',
    'HGNC' : 'hgnc_id', 
    'UniProtKB-ID' : 'uniprot_id', 
    'GeneID' : 'entrez_gene_id',

    # Entrez Terms
    'entrez_gene_id' : 'entrez_gene_id',
    'ensembl_gene_id' : 'ensembl_gene_id',
}

## GeneMANIA Settings

# Path to GeneMANIA ID Mapping File
PATH_TO_GENEMANIA_ID_MAP = PATH_TO_NOMENCLATURE + 'genemania/'
GENEMANIA_ID_MAP_FILENAME = 'genemania.human.id.map'

# List of Desired Fields
GENEMANIA_ID_MAP_FIELDS = [
    'Ensembl Gene ID',
    'Ensembl Protein ID',
    'Entrez Gene ID',
    'Gene Name',
    'RefSeq Protein ID',
    'RefSeq mRNA ID',
    'Synonym',
    'Uniprot ID',
    'Primary Gene ID',
    ] 

## BioGRID Settings

# Path to BioGRID ID Mapping File
PATH_TO_BIOGRID_ID_MAP = PATH_TO_NOMENCLATURE + 'biogrid/'
RAW_BIOGRID_ID_MAP_FILENAME = 'raw.biogrid.all.id.map'
BIOGRID_ID_MAP_FILENAME = 'biogrid.human.id.map'

# List of Desired Fields
BIOGRID_ID_MAP_FIELDS = [
    'ENTREZ_GENE',
    'SYNONYM',
    'REFSEQ_PROTEIN_ACCESSION',
    'ENSEMBL',
    'REFSEQ_DNA_ACCESSION',
    'OFFICIAL_SYMBOL',
    'ENTREZ_GENE',
    'HGNC',
]

## UniProt Settings

# Path to UniProt ID Mapping File
PATH_TO_UNIPROT_ID_MAP = PATH_TO_NOMENCLATURE + 'uniprot/'
UNIPROT_ID_MAP_FILENAME = 'uniprot.human.id.map'

# List of Desired Fields
UNIPROT_ID_MAP_FIELDS = [
    'RefSeq',
    'Ensembl_PRO',
    'Ensembl',
    'UniGene',
    'HGNC',
    'UniProtKB-ID',
    'GeneID',
]

## Entrez Settings

# Path to Entrez ID Mapping File
PATH_TO_ENTREZ_ID_MAP = PATH_TO_NOMENCLATURE + 'ncbi/'
RAW_ENTREZ_ID_MAP_FILENAME = 'raw.entrez.all.id.map'
ENTREZ_ID_MAP_FILENAME = 'entrez.human.id.map'

# List of Desired Fields
ENTREZ_ID_MAP_FIELDS = [
    'entrez_gene_id',
    'ensembl_gene_id',
]


# - - - - - - - - - - GLOBAL ID MAP CREATION - - - - - - - - - - #


# getNCBIOfficialSymbol: Gets official symbol for <entrezID>
def getNCBIOfficialSymbol(entrezID):
    BASE_NCBI_URL = 'http://www.ncbi.nlm.nih.gov/gene/6634'
    
    return
    
# createGeneIDMap: Creates database <geneIDMap> which maps between various
# gene nomenclatures, including Ensembl, Entrez, and RefSeq.
def createGeneIDMap():
    processList = []

    # Create Individual ID Maps
    for source, function in GENE_ID_MAP_CREATION.iteritems():
        processList.append( Process(target = function) )

    # Start Processes
    for process in processList:
        process.start()

    # Join Processes
    for process in processList:
        process.join()
    
    processList = []

    # Merge ID Maps
    if len(GENE_ID_MAP_MERGING.keys()) > 1:
        for source, function in GENE_ID_MAP_MERGING.iteritems():
            processList.append( Process(target = function) )
            
        # Start Processes
        for process in processList:
            process.start()

        # Join Processes
        for process in processList:
            process.join()
    else:
        numRecords = MongoClient().db.entrezIDMap.count()
        parallelizeTaskByRecord(mergeEntrezWithGenemaniaIDMap, numRecords)

    return


# - - - - - - - - - - GENEMANIA ID MAP CREATION - - - - - - - - - - #


# createGenemaniaIDMap: Create <geneIDMap> from GeneMANIA mapping file.
def createGenemaniaIDMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap

    # Open Input File
    inFile = open(PATH_TO_GENEMANIA_ID_MAP + GENEMANIA_ID_MAP_FILENAME, 'r')

    # Skip First Line
    skipLine = True

    # Start with Empty Record
    record = {}
    currentGene = ''

    # Iterate through File
    for line in inFile:
        # Skip First Line
        if (skipLine):
            skipLine = False
            continue

        # Get Line Fields
        lineFields = parseTabSeparatedLine(line)

        # If this is a new gene...
        if (lineFields[0] != currentGene):
            # ...and record is not empty, commit to DB
            if record:
                geneIDMap.insert(record)
                record = {}
            currentGene = lineFields[0]
            record.update( { 'primary_gene_id' : [currentGene] } )
        
        # Add Line Data to Current Record
        fieldName = CANONICAL_GENE_FIELDS.get(lineFields[2])
        if (record.get(fieldName)):
            # Take Union of 2 Lists
            record[fieldName] = list(set(record[fieldName]) | 
                                     set([lineFields[1]]))
        else:
            record.update( { fieldName : [lineFields[1]] } )

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return

# - - - - - - - - - - BIOGRID ID MAP CREATION - - - - - - - - - - #


# parseRawBiogridIDMap: Converts original BioGRID ID map file to better
# format for database creation, ignoring copyright, non-human organisms, etc.
def parseRawBiogridIDMap():
    # Open Input File
    inFile = open(PATH_TO_BIOGRID_ID_MAP + RAW_BIOGRID_ID_MAP_FILENAME, 'r')

    # Open Output File
    outFile = open(PATH_TO_BIOGRID_ID_MAP + BIOGRID_ID_MAP_FILENAME, 'w')

    # Ignore Header
    for i in range(28):
        inFile.readline()

    # Get Column Names
    headerLine = inFile.readline()
    outFile.write(headerLine)

    # Write Human Lines to Output
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        organism = lineFields[3]
        if (organism == DESIRED_ORGANISM):
            outFile.write(line)

    # Close Files
    inFile.close()
    outFile.close()

    return

# createBiogridIDMap: Create <biogridIDMap> from ID mapping file.
def createBiogridIDMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    coll = db.biogridIDMap

    # Open Input File
    inFile = open(PATH_TO_BIOGRID_ID_MAP + BIOGRID_ID_MAP_FILENAME, 'r')

    # Skip Column Info
    skipLine = True

    # Initialize Empty Record
    record = {}
    currentGene = ''

    # Iterate through File
    for line in inFile:
        # Skip Column Info
        if skipLine:
            skipLine = False
            continue
        
        lineFields = parseTabSeparatedLine(line)
        [biogridID, idValue, idType] = lineFields[0:3]

        # If this is a new gene...
        if (biogridID != currentGene):
            # ...and record is not empty, commit to DB
            if record:
                coll.insert(record)
                record = {}
            currentGene = biogridID
        
        # Add Line Data to Current Record
        if idType not in BIOGRID_ID_MAP_FIELDS:
            continue

        fieldName = CANONICAL_GENE_FIELDS.get(idType)
        if (record.get(fieldName)):
            # Take Union of 2 Lists
            record[fieldName] = list( set( record[fieldName] ) | 
                                      set( [idValue] ))
        else:
            record.update( { fieldName : [ idValue ] } )

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return


# - - - - - - - - - - UNIPROT ID MAP CREATION - - - - - - - - - - #


# createUniprotIDMap: Creates <uniprotIDMap> from UniProt ID mappings.
def createUniprotIDMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    coll = db.uniprotIDMap

    # Open Input File
    inFile = open(PATH_TO_UNIPROT_ID_MAP + UNIPROT_ID_MAP_FILENAME, 'r')
    
    # Initialize Empty Record
    record = {}
    currentGene = ''

    # Iterate through File
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        [uniprotID, idType, idValue] = lineFields[0:3]

        # If this is a new gene...
        if (uniprotID != currentGene):
            # ...and record is not empty, commit to DB
            if record:
                coll.insert(record)
                record = {}
            currentGene = uniprotID
        
        # Add Line Data to Current Record
        if idType not in UNIPROT_ID_MAP_FIELDS:
            continue

        fieldName = CANONICAL_GENE_FIELDS.get(idType)
        
        if fieldName == 'refseq_protein_id':
            idValue = idValue.split('.')[0]

        if fieldName == 'hgnc_id':
            idValue = idValue.split(':')[1]

        if (record.get(fieldName)):
            # Take union of 2 lists
            record[fieldName] = list( set( record[fieldName] ) | 
                                      set( [idValue] ))
        else:
            record.update( { fieldName : [ idValue ] } )

    # Close File and DB Connection
    inFile.close()
    cli.close()

    return


# - - - - - - - - - - ENTREZ ID MAP CREATION - - - - - - - - - - #


# parseRawEntrezIDMap: Converts original Entrez ID map file to better format
# for database creation, ignoring non-human organisms.
def parseRawEntrezIDMap():
    # Open Input File
    inFile = open(PATH_TO_ENTREZ_ID_MAP + RAW_ENTREZ_ID_MAP_FILENAME, 'r')
    
    # Open Output File
    outFile = open(PATH_TO_ENTREZ_ID_MAP + ENTREZ_ID_MAP_FILENAME, 'w')

    # Iterate through File
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

# createEntrezIDMap: Creates <entrezIDMap> from ID mapping file.
def createEntrezIDMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    entrezIDMap = db.entrezIDMap

    # Convert Raw Mapping File
    parseRawEntrezIDMap()

    # Open Input File
    inFile = open(PATH_TO_ENTREZ_ID_MAP + ENTREZ_ID_MAP_FILENAME, 'r')

    # Iterate through Mapping File
    for line in inFile:
        # Skip Comments
        if line[0] == '#':
            continue

        lineFields = parseTabSeparatedLine(line)

        [ entrezID, ensemblID ] = lineFields[1:3]

        record = { 
            'entrez_gene_id' : entrezID,
            'ensembl_gene_id' : ensemblID,
            }

        entrezIDMap.insert(record)
        
    # Close File and DB Connection
    inFile.close()
    cli.close()

    return


# - - - - - - - - - - ID MAP MERGING - - - - - - - - - - #


# mergeBiogridWithGenemaniaIDMap: Merges <biogridIDMap> with <geneIDMap>.
def mergeBiogridWithGenemaniaIDMap(skipVal = 0, limitVal = 10000):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap
    biogridIDMap = db.biogridIDMap

    # Iterate through all Records of <biogridIDMap>
    for biogridRecord in biogridIDMap.find(timeout = False, skip = skipVal, 
                                           limit = limitVal):
        
        # Get Records that Match <biogridRecord>
        query = constructQueryForGeneRecord(biogridRecord)
        geneResultSet = geneIDMap.find(query)

        if geneResultSet.count() == 1:
            # Merge BioGRID and GeneMANIA Records
            mergedRecord = mergeMongoRecords(geneResultSet[0], biogridRecord)
        elif geneResultSet.count() == 0:
            # Add BioGRID Record to <geneIDMap>
            del biogridRecord['_id']
            entrezID = biogridRecord.get('entrez_gene_id')
            biogridRecord.update( { 'primary_gene_id' : entrezID } )
            geneIDMap.insert(biogridRecord)

    # Close DB Connection
    cli.close()

    return

# mergeUniprotWithGenemaniaIDMap: Merges <uniprotIDMap> with <geneIDMap>.
def mergeUniprotWithGenemaniaIDMap(skipVal = 0, limitVal = 10000):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap
    uniprotIDMap = db.uniprotIDMap
    
    # Iterate through all Records of <uniprotIDMap>
    for uniprotRecord in uniprotIDMap.find(timeout=False, skip = skipVal, 
                                           limit = limitVal):
        
        # Get Records that Match <uniprotRecord>
        query = constructQueryFromGeneRecord(uniprotRecord)
        geneResultSet = geneIDMap.find(query)

        if geneResultSet.count() == 1:
            # Merge UniProt and GeneMANIA Records
            mergedRecord = mergeMongoRecords(geneResultSet[0], uniprotRecord)
        elif geneResultSet.count() == 0:
            # Add UniProt Record to <geneIDMap>
            del uniprotRecord['_id']
            uniprot = uniprotRecord.get('uniprot_id')
            uniprotRecord.update( { 'primary_gene_id' : uniprot } )
            geneIDMap.insert(uniprotRecord)

    # Close DB Connection
    cli.close()

    return

# mergeEntrezWithGenemaniaIDMap: Merges <entrezIDMap> with <geneIDMap>.
def mergeEntrezWithGenemaniaIDMap(skipVal = 0, limitVal = 20000):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap
    entrezIDMap = db.entrezIDMap

    # Iterate through <entrezIDMap>
    for entrezRecord in entrezIDMap.find(skip = skipVal, limit = limitVal):
        entrezID = entrezRecord.get('entrez_gene_id')
        ensemblID = entrezRecord.get('ensembl_gene_id')

        geneRecord = geneIDMap.find_one( { 'entrez_gene_id' : entrezID } )
        
        if geneRecord:
            # Add Ensembl Mapping
            if geneRecord.get('ensembl_gene_id'):
                if ensemblID in geneRecord['ensembl_gene_id']:
                    continue
                else:
                    geneRecord['ensembl_gene_id'].append(ensemblID)
        else:
            # Add Entrez Mapping
            geneRecord = geneIDMap.find_one( { 'ensembl_gene_id' : ensemblID } )
            if geneRecord:
                if geneRecord.get('entrez_gene_id'):
                    if entrezID in geneRecord['entrez_gene_id']:
                        continue
                    else:
                        geneRecord['entrez_gene_id'].append(entrezID)
            else:
                continue
            
        geneIDMap.save(geneRecord)

    # Close DB Connection
    cli.close()

    return


# - - - - - - - - - - FUNCTION MAPS - - - - - - - - - - #


GENE_ID_MAP_CREATION = { 
    'genemania' : createGenemaniaIDMap,
    # 'biogrid' : createBiogridIDMap,
    # 'uniprot' : createUniprotIDMap,
    'entrez' : createEntrezIDMap,
    }

GENE_ID_MAP_MERGING = { 
    #'biogrid' : mergeBiogridWithGenemaniaIDMap,
    # 'uniprot' : mergeUniprotWithGenemaniaIDMap,
    'entrez' : mergeEntrezWithGenemaniaIDMap,
}


# - - - - - - - - - - ANALYSIS - - - - - - - - - - #


# analyzeGeneIDMapCoverage: Counts occurrences of each nomenclature per DB.
def analyzeIDMapCoverage():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    genemaniaIDMap = db.geneIDMap
    uniprotIDMap = db.uniprotIDMap
    biogridIDMap = db.biogridIDMap
    entrezIDMap = db.entrezIDMap

    outFile = open(PATH_TO_ANALYSIS + 'gene.id.map.coverage', 'w')

    dbDict = {
        # 'BioGRID' : biogridIDMap,
        # 'UniProt' : uniprotIDMap,
        'GeneMANIA' : genemaniaIDMap,
        'Entrez' : entrezIDMap
        }

    for dbName, dbHandle in dbDict.iteritems():
        if PRINT_PROGRESS:
            print dbName

        outFile.write('Analyzing %s ID Map...\n' % dbName)

        geneCounts = {}
        FIELD_LIST = []

        if dbName == 'BioGRID':
            FIELD_LIST = BIOGRID_ID_MAP_FIELDS
        elif dbName == 'UniProt':
            FIELD_LIST = UNIPROT_ID_MAP_FIELDS
        elif dbName == 'Entrez':
            FIELD_LIST = ENTREZ_ID_MAP_FIELDS
        elif dbName == 'GeneMANIA':
            FIELD_LIST = GENEMANIA_ID_MAP_FIELDS

        # Count Occurrences of Each Gene Type
        for field in FIELD_LIST:
            canonField = CANONICAL_GENE_FIELDS.get(field)
            geneCounts.update( { canonField : 0 } )

        for record in dbHandle.find(timeout=False):
            for key, value in record.iteritems():
                if key == '_id':
                    continue
                geneCounts[key] = geneCounts[key] + 1

        # Write Output
        for key, value in geneCounts.iteritems():
            outFile.write('\t%s: %d\n' % (key, value))
        outFile.write('\n')

    cli.close()

