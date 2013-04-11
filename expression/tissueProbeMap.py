#! /usr/bin/python

# tissueProbeMap.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: March 24, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.databases import *
from common.statistics import *


# - - - - - - - - - - TISSUE:PROBE MAP CREATION - - - - - - - - - - #


# createTissueProbeMap: Creates <tissueProbeMap> from <sampleProbeMap>.
def createTissueProbeMap(minimum = 0, maximum = 80):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    sampleProbeMap = db.sampleProbeMap
    tissueProbeMap = db.tissueProbeMap

    # Iterate through all Tissues
    for tissue in FUNCTIONAL_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue

        # Initialize Record
        record = {
            'tissue' : tissue,
            'geo_id' : [],
            'probes' : {},
            }

        # Entry in record['probes'] = {
        #                             'probe_id' : [],
        #                             }
        query = constructQueryForTissue(tissue)
        sampleRecordList = sampleProbeMap.find(query)

        geoIDs = record['geo_id']
        probes = record['probes']

        # Iterate through Matching Samples
        for sampleRecord in sampleRecordList:
            # Add GEO ID
            geoIDs.append( sampleRecord['geo_id'] )
            # Add Probes from Sample
            for probeID, expression in sampleRecord['probes'].iteritems():
                if probeID not in probes:
                    probes[probeID] = list()
                probes[probeID].append( expression )
                
        tissueProbeMap.save( record )

    # Close DB Connection
    cli.close()
    
    return

# normalizeTissueProbeMap: Normalizes <tissueProbeMap> by computing median.
def normalizeTissueProbeMap(minimum = 0, maximum = 80):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    tissueProbeMap = db.tissueProbeMap
    normTissueProbeMap = db.normTissueProbeMap

    # Iterate through all Tissues
    for tissueRecord in tissueProbeMap.find()[minimum:maximum]:
        tissue = tissueRecord.get('tissue')
        if PRINT_PROGRESS:
            print tissue

        # Initialize Record
        record = {
            'geo_id' : tissueRecord.get('geo_id'),
            '_id' : tissueRecord.get('_id'),
            'tissue' : tissue,
            'probes' : {},
            }

        # Iterate through Probes
        for probeID, expressionList in tissueRecord['probes'].iteritems():
            expressionVals = [ float(s) for s in expressionList ]
            expressionValue = median(expressionVals)
            record['probes'].update( { probeID : expressionValue } )

        normTissueProbeMap.save(record)

    # Close DB Connection
    cli.close()

    return


# - - - - - - - - - - ANALYSIS - - - - - - - - - - #


# analyzeTissueProbeMap: Computes Statistics on Probe Expression Values.
def analyzeTissueProbeMap(verbose=False):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    tissueProbeMap = db.tissueProbeMap

    # Open Output File
    outFile = open(PATH_TO_EXPRESSION_ANALYSIS + 'tissue.probe.map.analysis', 'w')
    outFile.write('Analyzing Tissue:Probe Map...\n\n')

    # Iterate through DB
    for tissueRecord in tissueProbeMap.find():
        outFile.write('Tissue: %s\n' % tissueRecord['tissue'])
        outFile.write('\tSamples: %s\n' % str(tissueRecord['geo_id']))
        outFile.write('\tProbes:\n')
        
        # Track Expressed, Unexpressed, and Ambiguous Probes
        expressedProbes = []
        unexpressedProbes = []
        ambiguousProbes = {}
        
        for probeID, expressionList in tissueRecord['probes'].iteritems():
            expressionValues = [float(s) for s in expressionList] 

            # Count # Expressed Probes
            numExpressed = 0
            for value in expressionValues:
                if value > ARRAY_EXPRESSION_THRESHOLD:
                    numExpressed = numExpressed + 1

            # If All Expressed, Add to <expressedProbes>
            if numExpressed == len(expressionValues):
                expressedProbes.append(probeID)
            # If None Expressed, add to <unexpressedProbes>
            elif numExpressed == 0:
                unexpressedProbes.append(probeID)
            else:
                ambiguousProbes.update( { probeID : expressionValues } )

        # Write Summary to File
        outFile.write('\t\tExpressed Probes: %s\n' % len(expressedProbes))
        outFile.write('\t\tUnexpressed Probes: %s\n' % len(unexpressedProbes))
        outFile.write('\t\tAmbiguous Probes: %s\n' % len(ambiguousProbes))

        for probeID, expressionList in ambiguousProbes.iteritems():
            outFile.write('\t\t\t%s\t%s\t' % (probeID ,str(expressionList)))
            medi = median(expressionList)
            average, sem = mean(expressionList)
            dev = stdDev(expressionList)
            vote = poll(expressionList, ARRAY_EXPRESSION_THRESHOLD)
            outFile.write('%f\t%f\t%f\t%f\t%f\n' % (medi, average, sem, dev, vote))
                    
        outFile.write('\n\n')

    # Close File and DB Connection
    outFile.close()
    cli.close()
    
    return

