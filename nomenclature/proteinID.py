#! /usr/bin/python

# Author: Santhosh Balasubramanian
# Created: March 30, 2013
# Last Modified: March 30, 2013


# Library Imports
from pymongo import *


# getRefseqProteinForEnsemblGene: Returns Refseq Protein ID for Ensembl Gene ID.
def getRefseqProteinForEnsemblGene(ensemblID):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap
    uniprotDB = db.uniprotIDMap

    # Try Uniprot DB
    proteinRecord = uniprotDB.find_one( { 'ensembl_gene_id' : ensemblID } )

    if proteinRecord:
        refseqList = proteinRecord.get('refseq_protein_id')
        if refseqList:
            refseqID = refseqList[0]
        else:
            refseqID = None
    else:
        # Get Entrez ID
        geneRecord = geneIDMap.find_one( { 'ensembl_gene_id' : ensemblID } )
        entrezList = geneRecord.get('entrez_gene_id')
        if entrezList:
            entrezID = entrezList[0]
            refseqID = getRefseqProteinForEntrezGene(entrezID)
        else:
            refseqID = None

    cli.close()
    return refseqID
    
# getRefseqProteinForEntrezGene: Returns Refseq Protein ID for Entrez Gene ID.
def getRefseqProteinForEntrezGene(entrezID):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    uniprotDB = db.uniprotIDMap

    proteinRecord = uniprotDB.find_one( { 'entrez_gene_id' : entrezID } )
    if proteinRecord:
        refseqList = proteinRecord.get('refseq_protein_id')
        if refseqList:
            refseqID = refseqList[0]
        else:
            refseqID = None
    else:
        refseqID = None

    cli.close()
    return refseqID
