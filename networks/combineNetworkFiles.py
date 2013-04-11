#! /usr/bin/python/

# combineNetworkFiles.py
# Author: Santhosh Balasubramanian
# Date: January 27, 2013

import os

FILE_PREFIXES = {
    'ppi' : 'Physical_interactions',
    'genetic' : 'Genetic_interactions',
    'predicted' : 'Predicted',
    'coexpression' : 'Co-expression',
    'colocalization' : 'Co-localization',
    'pathway' : 'Pathway',
    'protein_domains' : 'Shared_protein_domains',
    }

PATH_TO_INPUT = '/home/santhosh/Dropbox/Thesis/data/networks/genemania/Homo_sapiens/'
PATH_TO_OUTPUT = '/home/santhosh/Dropbox/Thesis/data/networks/genemania/Homo_sapiens.COMBINED/'
OUTPUT_BASE_FILENAME = 'combined.%s.network'

def combineGenemaniaNetworks():
    
    # List all network files in directory
    inputDirFiles = os.listdir(PATH_TO_INPUT)
    
    # Iterate through each interaction type for combination file
    for interactionType, prefix in FILE_PREFIXES.iteritems():
        print interactionType
        interactionFiles = []

        # Get all network files for this interaction type
        for filename in inputDirFiles:
            if prefix in filename:
                interactionFiles.append(filename)

        # Add each interactionFile to outputfile for interaction type
        outputFilename = OUTPUT_BASE_FILENAME % interactionType
        outputFile = open(PATH_TO_OUTPUT + outputFilename, 'w')
        
        interactionSet = set()
        
        outputFile.write('Gene_A\tGene_B\tWeight\n')
        
        # Iterate through each interaction file
        for filename in interactionFiles:
            f = open(PATH_TO_INPUT + filename, 'r')
            # Skip header line
            f.readline()
            
            # Add lines of f to outputFile
            for line in f:
                lineTabSplit = line.split('\t')
                lineTabSplit = [s.strip() for s in lineTabSplit]
                geneA = lineTabSplit[0]
                geneB = lineTabSplit[1]
                # Don't write same edge to outputfile twice
                if (geneA, geneB) in interactionSet:
                    continue
                if (geneB, geneA) in interactionSet:
                    continue
                interactionSet.add((geneA, geneB))
                outputFile.write(line)

            f.close()

        outputFile.close()

