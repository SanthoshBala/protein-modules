#! /usr/bin/python

# ncbiBlast.py
# Author: Santhosh Balasubramanian
# Created: March 30, 2013
# Last Modified: March 30, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.strings import *
from common.statistics import *

# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

NCBI_REFSEQ_BASE_STRING = 'ref|%s|'
PATH_TO_SIMILAR_PROTEIN_ANALYSIS = PATH_TO_MODULE_SIMILARITY + 'exchange/'
PATH_TO_EXCHANGE_SIMILARITY = PATH_TO_MODULE_SIMILARITY + 'exchange/'


# - - - - - - - - - - FUNCTIONS - - - - - - - - - - #


# createNCBIBlastInput: Creates all input files for NCBI Blast. <numComponents>
# is number of connected components in debruijn graph. <numShards> is number
# of random divisions of protein lists.
def createNCBIBlastInput(numComponents, numShards):
    for i in range(numComponents):
        inFileName = 'debruijn.exchange.proteins.%d' % (i + 1)
        convertRefSeqProteinListToNCBIFormat(PATH_TO_EXCHANGE_SIMILARITY, \
                                             inFileName)

    for i in range(numShards):
        inFileName = 'debruijn.exchange.proteins.shard.%d' % (i + 1)
        convertRefSeqProteinListToNCBIFormat(PATH_TO_EXCHANGE_SIMILARITY, \
                                                 inFileName)
    
    return


# convertRefSeqProteinListToNCBIFormat: Writes a line-separated list of 
# Refseq Protein IDs in NCBI format for use with Blast.
def convertRefSeqProteinListToNCBIFormat(inFilePath, inFileName):
    # Open Input File
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    outFileName = inFileName + '.ncbi'
    outFilePath = inFilePath
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through <inFile>
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        proteinID = lineFields[0]

        ncbiString = NCBI_REFSEQ_BASE_STRING % proteinID

        outFile.write('%s\n' % ncbiString)

    # Close Files
    inFile.close()
    outFile.close()

    return

# parseNcbiBlastOutput: Parses NCBI blast output, writing results to DB.
def parseNcbiBlastOutput(inFilePath, inFileName):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    seqDB = db.proteinAlignment

    # Open Input File
    inFile = open(inFilePath + inFileName, 'r')

    currentRawQueryID = ''
    currentRawSubjectID = ''
    cleanQueryID = ''
    cleanSubjectID = ''

    firstRecord = True

    # Iterate through File
    for line in inFile:
        if line[0] == '\n':
            continue
        if line[0] == '#':
            continue
        
        lineFields = parseTabSeparatedLine(line)

        # Check Query Protein
        rawQueryID = lineFields[0]
        if rawQueryID != currentRawQueryID:
            currentRawQueryID = rawQueryID
            cleanQueryID = rawQueryID.split('|')[-2].split('.')[0]
            cleanQueryID = cleanQueryID.split('.')[0]

            # Add to DB
            if firstRecord:
                firstRecord = False
            else:
                seqDB.save(record)
            
            record = { 
                'refseq_protein_id' : cleanQueryID,
                'alignment' : [],
                }
            
        
        # Check Subject Protein
        rawSubjectID = lineFields[1]
        if rawSubjectID == currentRawSubjectID:
            continue
        else:
            currentRawSubjectID = rawSubjectID
            cleanSubjectID = rawSubjectID.split(';')[0].split('|')[-2]
            cleanSubjectID = cleanSubjectID.split('.')[0]

            # Get E Value
            eValue = lineFields[-2]

            # Add to Record
            record['alignment'].append( ( cleanSubjectID, eValue ) )

    seqDB.save(record)
    # Close File and DB Connection
    inFile.close()
    cli.close()

    return
    
# parseAllNcbiBlastOutput: Parses all NCBI blast output.
def parseAllNcbiBlastOutput():
    BASE_FILENAME = 'debruijn.exchange.proteins.%s.blast'
    inFilePath = PATH_TO_SIMILAR_PROTEIN_ANALYSIS

    for i in range(1,21):
        inFileName = BASE_FILENAME % str(i)
        parseNcbiBlastOutput(inFilePath, inFileName)

# getSequenceHomologyHistogram: Writes histogram showing how many protein
# pairs are homologous proteins.
def getSequenceHomologyHistogram():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    seqDB = db.proteinAlignment

    # Open Input File
    inFileName = 'debruijn.exchange.protein.pairs'
    inFilePath = PATH_TO_SIMILAR_PROTEIN_ANALYSIS
    inFile = open(inFilePath + inFileName, 'r')

    # Open Output File
    outFileName = 'sequence.homology.histogram'
    outFilePath = PATH_TO_SIMILAR_PROTEIN_ANALYSIS
    outFile = open(outFilePath + outFileName, 'w')

    eValueList = []

    # Iterate through Pairs
    for line in inFile:
        lineFields = parseTabSeparatedLine(line)
        
        proteinA = lineFields[0]
        proteinB = lineFields[1]
        if PRINT_PROGRESS:
            print proteinA

        recordA = seqDB.find_one( { 'refseq_protein_id' : proteinA } )
        if not recordA:
            continue

        # Iterate through A's Alignments
        alignments = recordA.get('alignment')
        for alignment in alignments:
            if alignment[0] == proteinB:
                eValueList.append(float(alignment[1]))

    # Write to File
    for value in eValueList:
        outFile.write('%f\n' % value)

    # Close Files and DB Connection
    inFile.close()
    outFile.close()
    cli.close()
    
    return
        
    
