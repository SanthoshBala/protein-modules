#! /usr/bin/python

# probeGeneMap.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: March 24, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.strings import *
from common.databases import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


GENE_NOMENCLATURES = [ 
        'primary_gene_id', 
        'ensembl_gene_id', 
        'entrez_gene_id',
        'gene_symbol',
        'unigene_id',
        ]

## Affymetrix Settings

PATH_TO_AFFY_ANNOTATION = PATH_TO_EXPRESSION + 'gse1133/affymetrix/'
AFFY_ANNOTATION_FILENAME = 'u133a.probe.annotation.csv'

# Maps Affy field names to canonical values.
CANONICAL_AFFY_FIELDS = {
    'Probe Set ID' : 'probe_set_id',
    'GeneChip Array' : 'genechip_array',
    'UniGene ID' : 'unigene_id',
    'Gene Title' : 'gene_title',
    'Gene Symbol': 'gene_symbol',
    'Chromosomal Location' : 'chromosomal_location',
    'Unigene Cluster Type' : 'unigene_cluster_type',
    'Ensembl' : 'ensembl_gene',
    'Entrez Gene' : 'entrez_gene',
    'RefSeq Protein ID' : 'refseq_protein_id',
    'RefSeq Transcript ID' : 'refseq_transcript_id',
    'Transcript ID(Array Design)' : 'unigene_id',
    'Target Description' : 'target_description'
}

## GNF Settings

PATH_TO_GNF_ANNOTATION = PATH_TO_EXPRESSION + 'gse1133/biogps/'
GNF_ANNOTATION_FILENAME = 'gnf1h.probe.annotation.tsv'

# Maps GNF field names to canonical values.
CANONICAL_GNF_FIELDS = {
    'ProbesetID' : 'probe_set_id',
    'RefSeq' : 'refseq_transcript_id',
    'UniGene' : 'unigene_id',
    'EntrezGene' : 'entrez_gene',
    'Symbol' : 'gene_symbol',
    'Description' : 'gene_title',
    'Ensembl_representative' : 'ensembl_transcript',
}


# - - - - - - - - - - GLOBAL PROBE:GENE MAP CREATION - - - - - - - - - - #


# createProbeGeneMap: Creates probe:gene map from Affymetrix and GNF arrays.
def createProbeGeneMap():

    createAffyProbeGeneMap()
    createGnfProbeGeneMap()

    return


# - - - - - - - - - - AFFY PROBE:GENE MAP CREATION - - - - - - - - - - - #


# createAffyProbeGeneMap: Creates probe:gene map from Affymetrix array.
def createAffyProbeGeneMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    probeGeneMap = db.probeGeneMap

    # Open Input File
    inFile = open(PATH_TO_AFFY_ANNOTATION + AFFY_ANNOTATION_FILENAME, 'r')

    # Get Indices for Desired Fields
    desiredFieldIndices = {}

    # Read Header from First Line
    readHeader = True

    # Iterate through File
    for line in inFile:
        # Ignore Comments
        if line[0] == '#':
            continue

        lineFields = parseCommaSeparatedLine(line)
        
        # Parse Header
        if readHeader:
            for i in range(len(lineFields)):
                if lineFields[i] in CANONICAL_AFFY_FIELDS.keys():
                    desiredFieldIndices.update( { lineFields[i] : i } )
                else:
                    desiredFieldIndices.update( { lineFields[i] : -1 } )
            readHeader = False
            continue

        # Add Probe Info to DB
        record = {}
        for field in CANONICAL_AFFY_FIELDS.keys():
            index = desiredFieldIndices.get(field)
            fieldValue = lineFields[index]
            retVal = sanitizeAffyLineItem(field, fieldValue, record)

        probeGeneMap.insert(record)

    # Close File and DB Connection
    inFile.close()
    cli.close()
    
    return            
        
# sanitizeAffyLineItem: Cleans data from annotation file for DB.
def sanitizeAffyLineItem(fieldName, lineItem, recordDict):
    items = lineItem.split("///")
    items = [s.strip() for s in items]
    
    # Convert to Lowercase and Underscore
    canonicalField = CANONICAL_AFFY_FIELDS.get(fieldName)

    # Handle UniGene ID Field
    if canonicalField == 'unigene_id':
        if 'Hs.' not in items[0]:
            recordDict.update( { canonicalField : [] } )
            return 1

    # Handle Description Field
    if canonicalField == 'target_description':
        items = lineItem.split()
        items = [s.strip() for s in items]
        value = ''
        for item in items:
            # Parse out UniGene
            if '/ug' in item.lower():
                value = item[4:]
                if recordDict.get('unigene_id'):
                    recordDict['unigene_id'].append(value)
                else:
                    recordDict.update( { 'unigene_id' : [value] } )
            if '/gen' in item.lower():
                value = item[5:]
                if recordDict.get('gene_symbol'):
                    recordDict['gene_symbol'].append(value)
                else:
                    recordDict.update( { 'gene_symbol' : [value] } )
        return 1

    if (len(items) == 1) and (items[0] == '---'):
        items = []

    if recordDict.get(canonicalField):
        recordDict[canonicalField] = recordDict[canonicalField] + items
    else:
        recordDict.update( { canonicalField : items } )

    return 1


# - - - - - - - - - - - GNF PROBE:GENE MAP CREATION - - - - - - - - - - #


# createGnfProbeGeneMap: Create probe:gene map from GNF array.
def createGnfProbeGeneMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    probeGeneMap = db.probeGeneMap

    # Open Input File
    inFile = open(PATH_TO_GNF_ANNOTATION + GNF_ANNOTATION_FILENAME, 'r')

    # Get Indices for Desired Fields
    desiredFieldIndices = {}

    # Read Header from First Line
    readHeader = True

    # Iterate through File
    for line in inFile:
        # Ignore Comments
        if (line[0] == '#') or (line[0] == '!') or (line[0] == '^'):
            continue

        lineFields = parseTabSeparatedLine(line)

        # Read Header
        if readHeader:
            for i in range( len(lineFields) ):
                # Get Desired Field Indices
                if lineFields[i] in CANONICAL_GNF_FIELDS.keys():
                    desiredFieldIndices.update( { lineFields[i] : i } )
                else:
                    desiredFieldIndices.update( { lineFields[i] : -1 } )
                    
            readHeader = False
            continue
        
        # Add Probe Data to DB
        record = {}
        for field in CANONICAL_GNF_FIELDS.keys():
            index = desiredFieldIndices.get(field)
            fieldValue = lineFields[index]
            
            if field == 'Description':
                # If "Gene:" in Description, Parse out ensembl_gene
                if 'Gene:' in field:
                    descriptionFields = lineItem.split()
                    geneField = descriptionFields[0]
                    ensembleGene = geneField.split(':')[-1]
                    record.update( { 'ensembl_gene' : [ ensemblGene ] } )
                    continue

            canonicalField = CANONICAL_GNF_FIELDS.get(field)

            if fieldValue == '':
                value = []
            else:
                value = [ fieldValue ]

            record.update( { canonicalField : value } )
            
        probeGeneMap.insert(record)

    return
        
        
# - - - - - - - - - - EXTENDED PROBE:GENE MAP CREATION - - - - - - - - - - #


# createExtProbeGeneMap: Creates <extProbeGeneMap> from <probeGeneMap>, <geneIDMap>.
def createExtProbeGeneMap(skipVal = 0, limitVal = 10000):
    # Open db connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap
    probeGeneMap = db.probeGeneMap
    extProbeGeneMap = db.extProbeGeneMap

    # Iterate through <geneIDMap>
    for gene in geneIDMap.find(timeout=False,skip=skipVal,limit=limitVal):
        # Query for Probes Matching <gene>
        query = constructQueryForGeneRecord(gene)
        probeResultSet = probeGeneMap.find(query)

        # Iterate through Probes in <probeResultSet>
        for probeRecord in probeResultSet:
            # Extend <probeRecord> and commit to <extProbeGeneMap>
            extendedProbeRecord = mergeMongoRecords(probeRecord, gene)
            extProbeGeneMap.save(extendedProbeRecord)
    
    # Close DB Connection
    cli.close()

    return

# mergePGMwithEPGM: Merge <probeGeneMap> with <extProbeGeneMap>.
def mergePGMwithEPGM(skipVal = 0, limitVal = 10000):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    probeGeneMap = db.probeGeneMap
    extProbeGeneMap = db.extProbeGeneMap

    for pgmRecord in probeGeneMap.find(skip=skipVal, limit=limitVal):
        # Look up in EPGM by _id
        pgmID = pgmRecord.get('_id')
        recordSet = extProbeGeneMap.find( { '_id' : pgmID } )

        if recordSet.count() == 0:
            extProbeGeneMap.save(pgmRecord)

    # Close DB Connection
    cli.close()
    
    return


# - - - - - - - - - - ANALYSIS - - - - - - - - - - #


# analyzeProbeGeneMapCoverage: Prints coverage of each nomenclature type.
def analyzeProbeGeneMapCoverage(extended=False):
    # Open DB connection
    cli = MongoClient()
    db = cli.db
    if extended:
        coll = db.extProbeGeneMap
    else:
        coll = db.probeGeneMap

    # Open Output File
    if extended:
        outFile = open(PATH_TO_EXPRESSION_ANALYSIS + 'extended.probe.gene.coverage', 'w')
        outFile.write('Analyzing Extended Probe:Gene Map...\n\n')
    if not extended:
        outFile = open(PATH_TO_EXPRESSION_ANALYSIS + 'probe.gene.coverage', 'w')
        outFile.write('Analyzing Probe:Gene Map...\n\n')


    # Track Unique Genes within each Nomenclature
    geneSets = {}
    for nomenclature in GENE_NOMENCLATURES:
        geneSets.update( { nomenclature : set() } ) 

    # Track Probes Lacking Gene ID
    probeSets = {}
    for nomenclature in GENE_NOMENCLATURES:
        probeSets.update( { nomenclature : set() } )

    # Iterate through Probes
    for probeRecord in coll.find():
        # Track if Gene Type Exists for <probeRecord>
        geneBools = {}
        for nomenclature in GENE_NOMENCLATURES:
            geneBools.update( { nomenclature : True } )

        for nomenclature in GENE_NOMENCLATURES:
            if nomenclature == 'primary_gene_id':
                if not probeRecord.get(nomenclature):
                    geneBools[nomenclature] = False
            else:
                if not probeRecord.get(nomenclature):
                    geneBools[nomenclature] = False

        # Update <probeSet> and <geneSet>
        for key, value in geneBools.iteritems():
            if value == False:
                probeSet = probeSets.get(key)
                probeSet.add(probeRecord['probe_set_id'][0])
            else:
                geneSet = geneSets.get(key)
                for gene in probeRecord.get(key):
                    geneSet.add(gene)
                    
        outFile.write('\n')

    # Store values about what percentage of probes are covered by a given
    # gene nomenclature. i.e. Ensembl can identify 80% of probes
    coverageDict = {}
    numProbes = float(coll.count())
    for nomenclature in GENE_NOMENCLATURES:
        probeSet = probeSets.get(nomenclature)
        coverage = 1.0 - (len(probeSet)/numProbes)
        coverageDict.update( { nomenclature : coverage } )

    outFile.write('Summary\n')
    outFile.write('\tNumber of Probes: %f\n' % numProbes)
    for nomenclature in GENE_NOMENCLATURES:
        geneSet = geneSets.get(nomenclature)
        geneCount = len(geneSet)
        coverage = coverageDict.get(nomenclature)
        outFile.write('\tNumber of %s: %d\n' % (nomenclature, geneCount))
        outFile.write('\t%s Coverage: %f\n' % (nomenclature, coverage))

    # Close DB Connection
    cli.close()

    return
