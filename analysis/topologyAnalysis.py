#! /usr/bin/python

# topologyAnalysis.py
# Author: Santhosh Balasubramanian
# Created: March 12, 2013
# Last Modified: April 17, 2013


# Global Imports
from settings import *

# Graph Imports
from graphs.graphIO import *
from graphs.graphUtil import *

# Module Imports
from modules.moduleTopology import *

# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #

PATH_TO_TOPOLOGY_ANALYSIS = PATH_TO_ANALYSIS + 'topology/'
PATH_TO_LCC_ANALYSIS = PATH_TO_TOPOLOGY_ANALYSIS + 'lcc/'
PATH_TO_DEGREE_ANALYSIS = PATH_TO_TOPOLOGY_ANALYSIS + 'degree/'
PATH_TO_PATH_LENGTH_ANALYSIS = PATH_TO_TOPOLOGY_ANALYSIS + 'path_length/'
PATH_TO_MODULE_TOPOLOGY_ANALYSIS = PATH_TO_TOPOLOGY_ANALYSIS + 'module_topology/'

# - - - - - - - - - - MEAN PATH LENGTH - - - - - - - - - - #

# getNetworkMeanPathLengths
def getNetworkMeanPathLengths(minimum = 0, maximum = 85):

    # Iterate through All Tissues
    for tissue in AUGMENTED_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue

        # Open Output File
        outFileName = '%s.network.mean.path.length' % tissue
        outFilePath = PATH_TO_PATH_LENGTH_ANALYSIS
        outFile = open(outFilePath + outFileName, 'w')


        graph = getTissueSubgraph(tissue)
        distances = getPathLengthMatrix(graph)

        numVertices = graph.num_vertices()
        distSum = 0
        numPaths = 0
        # Sum over Shortest Paths
        for vertex in graph.vertices():
            distList = list(distances[vertex].a)
            for distance in distList:
                if PRINT_PROGRESS:
                    print '\t%f' % distance
                if distance > numVertices:
                    continue
                distSum = distSum + distance
                numPaths = numPaths + 1

        # Avoid Double Counting
        distSum = distSum / 2
        numPaths = numPaths / 2
        meanDist = distSum/float(numPaths)
        
        outFile.write('%s\t%f\n' % (tissue, meanDist))
        outFile.close()
    
    return

# - - - - - - - - - - DEGREE DISTRIBUTION - - - - - - - - - - #

# getNetworkDegreeDistributions: Gets degree distribution histogram.
def getNetworkDegreeDistributions(minimum = 0, maximum = 85):
    
    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue

        # Open Output File
        outFileName = '%s.network.degree.distribution' % tissue
        outFilePath = PATH_TO_DEGREE_ANALYSIS
        outFile = open(outFilePath + outFileName, 'w')

        graph = getTissueSubgraph(tissue)
        
        numVertices = graph.num_vertices()

        buckets, counts = getDegreeHistogram(graph)

        numCounts = len(counts)

        for i in range(numCounts):
            probability = counts[i]/float(numVertices)
            outFile.write('%f\t%f\n' % (buckets[i], probability))

        # Close File
        outFile.close()

    return


# - - - - - - - - - - LCC vs DEGREE - - - - - - - - - - #

def getNetworkLCCDegree(minimum = 0, maximum = 85):
    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue

        # Open Output File
        outFileName = '%s.network.degree.lcc' % tissue
        outFilePath = PATH_TO_LCC_ANALYSIS
        outFile = open(outFilePath + outFileName, 'w')

        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get LCC
        computeLocalClusteringCoefficients(graph)

        # Iterate through Graph
        lccProp = graph.vertex_properties['local_clustering_coefficient']
        for vertex in graph.vertices():
            degree = vertex.out_degree()
            lcc = lccProp[vertex]
            # Write to File
            outFile.write('%f\t%f\n' % (degree, lcc))

        # Close File
        outFile.close()

    return
                           
# - - - - - - - - - - MODULE TOPOLOGY ANNOTATION - - - - - - - - - - #

# getModuleTopologyAnnotation: For each tissue, gets the ontology annotation
# for all modules in module topology. Then prints a series of (x, y) coords
# for x = height in module topology and y = number of annotations
def getModuleTopologyAnnotation(minimum = 0, maximum = 85):
    
    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue
        
        # Open Output File
        functionFileName = '%s.module.topology.%s.annotation' % (tissue, 'function')
        processFileName = '%s.module.topology.%s.annotation' % (tissue, 'process')
        componentFileName = '%s.module.topology.%s.annotation' % (tissue, 'component')

        outFilePath = PATH_TO_MODULE_TOPOLOGY_ANALYSIS
        functionFile = open(outFilePath + functionFileName, 'w')
        processFile = open(outFilePath + processFileName, 'w')
        componentFile = open(outFilePath + componentFileName, 'w')

        # Get Module Topology
        graph = getTissueSubgraph(tissue)
        modTopology = getModuleTopology(graph)

        # Get Topology Annotation
        functions, processes, components = getTopologyAnnotation(modTopology)

        # Iterate through Vertices
        heightProp = modTopology.vertex_properties['height']
        geneIDProp = modTopology.vertex_properties['gene_id']
        for vertex in modTopology.vertices():
            geneID = geneIDProp[vertex]
            height = heightProp[vertex]
            
            geneFuncs = functions.get(geneID)
            if geneFuncs:
                numFunctions = len(geneFuncs)
            else:
                numFunctions = 0

            geneProcs = processes.get(geneID)
            if geneProcs:
                numProcesses = len(geneProcs)
            else:
                numProcesses = 0

            geneComps = components.get(geneID)
            if geneComps:
                numComponents = len(geneComps)
            else:
                numComponents = 0

            functionFile.write('%s\t%s\n' % (str(height), str(numFunctions)))
            processFile.write('%s\t%s\n' % (str(height), str(numProcesses)))
            componentFile.write('%s\t%s\n' % (str(height), str(numComponents)))
 
        functionFile.close()
        processFile.close()
        componentFile.close()

    return
