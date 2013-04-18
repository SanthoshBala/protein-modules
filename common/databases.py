#! /usr/bin/python

# database.py
# Author: Santhosh Balasubramanian
# Created: March 24, 2013
# Last Modified: March 24, 2013

# Python Imports
import json

# Library Imports
from pymongo import *

# Global Imports
from settings import *


# - - - - - - - - - - DATABASE CREATION - - - - - - - - - - #

# checkDatabaseDependencies: Checks that all database dependencies have been
# created for desiredDB. Generates dependencies if necessary.
def checkDatabaseDependencies(desiredDB):

    # INCOMPLETE

    return


# - - - - - - - - - - DATABASE DESTRUCTION - - - - - - - - - - #


# dropNomenclatureDBs: Drop all MongoDB databases that map gene nomenclatures.
def dropNomenclatureDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    
    # Drop Databases
    db.geneIDMap.drop()
    db.biogridIDMap.drop()
    db.uniprotIDMap.drop()
    db.entrezIDMap.drop()

    # Close DB Connection
    cli.close()
    return

# dropMicroarrayDBs: Drop all MongoDB databases that store intermediate
# results in processing of microarrays.
def dropMicroarrayDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db

    # Drop Databases
    dropProbeGeneMapDBs()


    db.tissueGeneMap.drop()
    db.normTissueGeneMap.drop()

    db.geneTissueMap.drop()
    db.normGeneTissueMap.drop()

    # Close DB Connection
    cli.close()

    return

# dropProbeGeneMapDBs: Drops <probeGeneMap>, <extProbeGeneMap>
def dropProbeGeneMapDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db

    # Drop Databases
    db.probeGeneMap.drop()
    db.extProbeGeneMap.drop()

    # Close DB Connection
    cli.close()

    return

# dropTissueProbeMapDBs: Drops <sampleProbeMap>, <tissueProbeMap>, and
# <normTissueProbeMap>.
def dropTissueProbeMapDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    
    # Drop Databases
    db.sampleProbeMap.drop()
    db.tissueProbeMap.drop()
    db.normTissueProbeMap.drop()

    # Close DB Connection
    cli.close()
    
    return

# dropGeneTissueMapDBs: Drops <tissueGeneMap>, <normTissueGeneMap>, 
# <geneTissueMap>, and <normGeneTissueMap>.
def dropGeneTissueMapDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db

    # Drop Databases
    db.tissueGeneMap.drop()
    db.normTissueGeneMap.drop()
    db.geneTissueMap.drop()
    db.normGeneTissueMap.drop()

    # Close DB Connection
    cli.close()
    
    return


# dropShuffleExpressionDBs: Drops <shuffleGeneTissueMap>,
# <shuffleNormGeneTissueMap>, and <shuffleNormTissueGeneMap>
def dropShuffleExpressionDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db

    # Drop Shuffle DBs
    db.shuffleGeneTissueMap.drop()
    db.shuffleNormGeneTissueMap.drop()
    db.shuffleNormTissueGeneMap.drop()

    # Close DB Connection
    cli.close()
    
    return
    
# dropModuleDBs: Drops <modules>.
def dropModuleDBs():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db

    # Drop Module DBs
    db.modules.drop()

    # Close DB Connection
    cli.close()
    
    return
    
# dropAllDBs: Drop all known MongoDB databases.
def dropAllDBs():

    # Drop Gene Nomenclature DBs
    dropNomenclatureDBs()

    # Drop Microarray DBs
    dropMicroarrayDBs()

    # Drop Probe:Gene Map DBs
    dropProbeGeneMapDBs()

    # Drop Tissue:Probe Map DBs
    dropTissueProbeMapDBs()

    # Drop Gene:Tissue Map DBS
    dropGeneTissueMapDBs()

    # Drop Shuffle Expression DBs
    dropShuffleExpressionDBs()

    return


# - - - - - - - - - - QUERYING - - - - - - - - - - #


# constructQueryForTissue: Constructs MongoDB Query.
def constructQueryForTissue(tissue):
    queryStr = '{"tissue" : "' + tissue + '"}'
    query = json.loads(queryStr)
    return query

# constructQueryForProbe
def constructQueryForProbe(probeID):
    queryStr = '{ "probe_set_id" : "%s" }' % probeID
    query = json.loads(queryStr)
    return query

# constructQueryForGeneRecord: Construct MongoDB Query that matches any of the
# nomenclature IDs in the record.
def constructQueryForGeneRecord(geneRecord):
    # Initialize Query
    recordList = '{"$or":['

    # Iterate through Record Keys
    for field in geneRecord.keys():
        # Ignore MongoDB ObjectID
        if field == '_id':
            continue

        # Construct List of Dictionaries of Values
        fieldList = '{"$or":['
        for value in geneRecord.get(field):
            fieldList = fieldList + '{"' + field + '":' + '"' + value + '"},'
        fieldList = fieldList[:-1] + ']},'

        recordList = recordList + fieldList

    recordList = recordList[:-1] + "]}"
    query = json.loads(recordList)

    return query

# constructQueryForGO: Construct MongoDB Query that matches the annotation.
def constructQueryForGO(goTerm, annotationType = 'function', geneList = None):
    # Initialize Query
    if geneList:
        orString = ','.join( [ '{ "entrez_gene_id" : "%s" }' % gene for gene in geneList ] )
        queryTemplateStr = '{ "$and" : [ { "%s" : [ "%s", "%s" ] }, { "$or" : [ %s ] } ] }'
        print queryTemplateStr
        queryStr = queryTemplateStr % (annotationType, goTerm[0], goTerm[1], 
                                       orString)
        print queryStr
    else:
        queryTemplateStr = '{ "%s" : [ "%s", "%s" ] }' 
        print queryTemplateStr
        queryStr = queryTemplateStr % (annotationType, goTerm[0], goTerm[1])
        print queryStr

    query = json.loads(queryStr)

    return query

# - - - - - - - - - - MISCELLANEOUS - - - - - - - - - - #


# mergeMongoRecords: Not a strictly generic dictionary union because
# assumes that value in all dictionaries are lists of values, and assumes
# the first argument keeps its ObjectID value.
def mergeMongoRecords(dictA, dictB):
    unionDict = {}

    # Add all values from dict A
    for key, value in dictA.iteritems():
        unionDict.update( { key : value } )

    for key, value in dictB.iteritems():
        # Ignore MongoDB's ObjectID
        if key == '_id':
            continue

        if unionDict.get(key):
            unionDict[key] = list( set(dictA[key]) | set(dictB[key]) )
        else:
            unionDict.update( { key : value } )

    return unionDict
