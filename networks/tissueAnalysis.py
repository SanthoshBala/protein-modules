#! /usr/bin/python

# Santhosh Balasubramanian
# March 12, 2013

from settings import *
from graphIO import *
from graphUtil import *
from moduleUtil import *

# - - - - - VERTICES - - - - - #

# getTissueVertexSetSizes
def getTissueVertexSetSizes():
    # Open Output File
    outFilename = 'tissue.ppi.num_vertices'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_TISSUE_LIST =  ['global'] + FUNCTIONAL_TISSUE_LIST
    
    # Iterate Through All Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        # Get Tissue Subgraph
        graph = getTissueSubgraph(tissue)

        # Get Vertex Set
        numVertices = getNumVertices(graph)

        # Write Size to File
        outFile.write('%s\t%d\n' % (tissue, numVertices))

    return

# getOrganSystemVertexSetSizes
def getOrganSystemVertexSetSizes():
    # Open Output File
    outFilename = 'organ_system.ppi.num_vertices'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    # Initialize System:Vertex Set Dict
    systemDict = {}
    for system in AUGMENTED_ORGAN_SYSTEM_LIST:
        systemDict.update( { system : set() } )

    # Iterate Through All Tissues
    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Get Tissue Subgraph
        graph = getTissueSubgraph(tissue)

        # Get Vertex Set
        vertexSet = getVertexSet(graph)

        # Get this Tissue's Organ System
        tissueSystem = TISSUE_SYSTEM_MAP.get(tissue)

        # Update both the Global and Organ System Vertex Sets
        systemDict[tissueSystem] = systemDict[tissueSystem].union(vertexSet)
        systemDict['global'] = systemDict['global'].union(vertexSet)

    # Write Sizes to File
    for system, vertexSet in systemDict.iteritems():
        numVertices = len(vertexSet)
        outFile.write('%s\t%d\n' % (system, numVertices))

    return

# - - - - - EDGES - - - - - #

# getTissueEdgeSetSizes
def getTissueEdgeSetSizes():
    # Open Output File
    outFilename = 'tissue.ppi.num_edges'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST
    
    # Iterate Through All Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        # Get Subgraph Filename
        filename = GENEMANIA_TISSUE_BASE_FILENAME % \
             (tissue, DESIRED_INTERACTION_TYPE)
        filepath = PATH_TO_TISSUE_SUBGRAPHS + filename

        # Get Num Edges
        numLines = getNumLinesInFile(filepath)
        # Subtract 1 to account for top line: "Gene_A Gene_B Weight"
        numEdges = numLines - 1

        # Write Size to File
        outFile.write('%s\t%d\n' % (tissue, numEdges))

    return

# getOrganSystemEdgeSetSizes
def getOrganSystemEdgeSetSizes():
    # Open Output File
    outFilename = 'organ_system.ppi.num_edges'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    # Initialize System:Vertex Set Dict
    systemDict = {}
    for system in AUGMENTED_ORGAN_SYSTEM_LIST:
        systemDict.update( { system : set() } )

    # Iterate Through All Tissues
    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Get Tissue Subgraph
        graph = getTissueSubgraph(tissue)

        # Get Edge Set
        edgeSet = getEdgeSet(graph)

        # Get this Tissue's Organ System
        tissueSystem = TISSUE_SYSTEM_MAP.get(tissue)

        # Update both the Global and Organ System Vertex Sets
        systemDict[tissueSystem] = systemDict[tissueSystem].union(edgeSet)
        systemDict['global'] = systemDict['global'].union(edgeSet)

    # Write Sizes to File
    for system, edgeSet in systemDict.iteritems():
        numEdges = len(edgeSet)
        outFile.write('%s\t%d\n' % (system, numEdges))

    return


# getEdgeTissueHistogram
def getEdgeTissueHistogram():
    # Open Output File
    outFilename = 'edge.tissue.histogram'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    histList = [0]*(len(AUGMENTED_TISSUE_LIST) + 1)

    # Initialize System:Vertex Set Dict
    systemDict = {}
    for system in AUGMENTED_TISSUE_LIST:
        systemDict.update( { system : set() } )

    # Iterate through all Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue

        # Get Tissue Subgraph
        graph = getTissueSubgraph(tissue)

        # Get Edge Set
        edgeSet = getEdgeSet(graph)

        # Update both the Global and Organ System Vertex Sets
        systemDict[tissue] = systemDict[tissue].union(edgeSet)
        systemDict['global'] = systemDict['global'].union(edgeSet)

    # Iterate through global
    for edge in systemDict['global']:
        count = 0
        for tissue in FUNCTIONAL_TISSUE_LIST:
            if edge in systemDict[tissue]:
                count = count + 1

        histList[count] = histList[count] + 1
        
    # Write to Output
    for i in range(len(histList)):
        outFile.write('%d\t%d\n' % (i, histList[i]))

    return

# getEdgeTissueHistogram
def getEdgeOrganSystemHistogram():
    # Open Output File
    outFilename = 'edge.organ_system.histogram'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST
    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    histList = [0]*(len(AUGMENTED_ORGAN_SYSTEM_LIST) + 1)

    # Initialize System:Vertex Set Dict
    systemDict = {}
    for system in AUGMENTED_ORGAN_SYSTEM_LIST:
        systemDict.update( { system : set() } )

    print systemDict

    # Iterate through all Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue

        # Get Tissue Subgraph
        graph = getTissueSubgraph(tissue)

        # Get Edge Set
        edgeSet = getEdgeSet(graph)

        # Get this Tissue's Organ System
        if tissue == 'global':
            tissueSystem = 'global'
        else:
            tissueSystem = TISSUE_SYSTEM_MAP.get(tissue)

        print tissue
        print tissueSystem
        # Update both the Global and Organ System Vertex Sets
        systemDict[tissueSystem] = systemDict[tissueSystem].union(edgeSet)
        systemDict['global'] = systemDict['global'].union(edgeSet)

    # Iterate through global
    for edge in systemDict['global']:
        count = 0
        for tissueSystem in CANONICAL_ORGAN_SYSTEM_LIST:
            if edge in systemDict[tissueSystem]:
                count = count + 1

        histList[count] = histList[count] + 1
        

    # Write to Output
    for i in range(len(histList)):
        outFile.write('%d\t%d\n' % (i, histList[i]))

    return


# - - - - - MODULES - - - - - #

# getTissueModuleSetSizes
def getTissueModuleSetSizes():
    # Open Output File
    outFilename = 'tissue.ppi.num_modules'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST
    
    # Iterate Through All Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue

        # Get Subgraph Filename
        filename = TISSUE_MODULE_IDS_BASE_FILENAME % tissue
        filepath = PATH_TO_TISSUE_MODULE_IDS + filename

        # Get Num Edges
        numModules = getNumLinesInFile(filepath)

        # Write Size to File
        outFile.write('%s\t%d\n' % (tissue, numModules))

    return

# getOrganSystemModuleSetSizes
def getOrganSystemModuleSetSizes():
    # Open Output File
    outFilename = 'organ_system.ppi.num_modules'
    outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    # Initialize System:Vertex Set Dict
    systemDict = {}
    for system in AUGMENTED_ORGAN_SYSTEM_LIST:
        systemDict.update( { system : set() } )

    # Iterate Through All Tissues
    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Get Module ID Set
        moduleSet = getTissueModuleIDSet(tissue)

        # Get this Tissue's Organ System
        tissueSystem = TISSUE_SYSTEM_MAP.get(tissue)

        # Update both the Global and Organ System Vertex Sets
        systemDict[tissueSystem] = systemDict[tissueSystem].union(moduleSet)
        systemDict['global'] = systemDict['global'].union(moduleSet)

    # Write Sizes to File
    for system, moduleSet in systemDict.iteritems():
        numModules = len(moduleSet)
        outFile.write('%s\t%d\n' % (system, numModules))

    return

# - - - - - DIAMETER - - - - - #

# getTissueSubgraphDiameters
def getTissueSubgraphDiameters():
    # Open Output File
    outFilename = 'tissue.ppi.diameter'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')
    
    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get Diameter
        diameter = getDiameter(graph)

        # Write to File
        outFile.write('%s\t%f\n' % (tissue, diameter))

    return

# getOrgansystemDiameters
def getOrganSystemDiameters():
    # Open Output File
    outFilename = 'organ_system.ppi.diameter'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST
    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    diameterDict = {}
    for organ in AUGMENTED_ORGAN_SYSTEM_LIST:
        diameterDict.update( { organ : [] } )

    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get Diameter
        diameter = getDiameter(graph)

        # Get Organ System
        if tissue == 'global':
            organSystem = 'global'
        else:
            organSystem = TISSUE_SYSTEM_MAP.get(tissue)

        # Write to Dict
        diameterDict[organSystem].append(diameter)

    for system, dList in diameterDict.iteritems():
        diameter = mean(dList)[0]
        outFile.write('%s\t%f\n' % (system, diameter))

    return
    

# - - - - - DENSITY - - - - - #

# getTissueSubgraphDensities
def getTissueSubgraphDensities():
    # Open Output File
    outFilename = 'tissue.ppi.density'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')
    
    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get Density
        density = getDensity(graph)

        # Write to File
        outFile.write('%s\t%f\n' % (tissue, density))

    return

# getOrganSystemDensities
def getOrganSystemDensities():
    # Open Output File
    outFilename = 'organ_system.ppi.density'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST
    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    densitiesDict = {}
    for organ in AUGMENTED_ORGAN_SYSTEM_LIST:
        densitiesDict.update( { organ : [] } )

    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get Density
        density = getDensity(graph)

        # Get Organ System
        if tissue == 'global':
            organSystem = 'global'
        else:
            organSystem = TISSUE_SYSTEM_MAP.get(tissue)

        # Write to Dict
        densitiesDict[organSystem].append(density)

    for system, dList in densitiesDict.iteritems():
        density = mean(dList)[0]
        outFile.write('%s\t%f\n' % (system, density))

    return

# - - - - - GLOBAL CLUSTERING COEFFICIENT - - - - - #

# getTissueGCCs()
def getTissueSubgraphGCCs():
    # Open Output File
    outFilename = 'tissue.ppi.gcc'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get GCC
        gcc = getGlobalClusteringCoefficient(graph)[0]

        # Write to File
        outFile.write('%s\t%f\n' % (tissue, gcc))

    return

# getOrganSystemGCCs
def getOrganSystemGCCs():
    # Open Output File
    outFilename = 'organ_system.ppi.gcc'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST
    AUGMENTED_ORGAN_SYSTEM_LIST = ['global'] + CANONICAL_ORGAN_SYSTEM_LIST

    gccDict = {}
    for organ in AUGMENTED_ORGAN_SYSTEM_LIST:
        gccDict.update( { organ : [] } )

    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get GCC
        gcc = getGlobalClusteringCoefficient(graph)[0]

        # Get Organ System
        if tissue == 'global':
            organSystem = 'global'
        else:
            organSystem = TISSUE_SYSTEM_MAP.get(tissue)

        # Write to Dict
        gccDict[organSystem].append(gcc)

    for system, gccList in gccDict.iteritems():
        gcc = mean(gccList)[0]
        outFile.write('%s\t%f\n' % (system, gcc))
    
    return

# - - - - - LOCAL CLUSTERING COEFFICIENT - - - - - #

def getTissueSubgraphLCCDegreePlots():
    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    for tissue in AUGMENTED_TISSUE_LIST:
        # Open Output File
        outFilename = '%s.ppi.degree_lcc' % tissue
        outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
        outFile = open(outFilepath, 'w')
        
        print tissue
        # Get Graph
        graph = getTissueSubgraph(tissue)

        # Get LCC
        computeLocalClusteringCoefficients(graph)

        # Iterate through graph
        lccProp = graph.vertex_properties['local_clustering_coefficient']
        for vertex in graph.vertices():
            degree = vertex.out_degree()
            lcc = lccProp[vertex]
            # Write to File
            outFile.write('%f\t%f\n' % (degree, lcc))

    return

# - - - - - MEAN PATH LENGTHS - - - - - #

def getTissueSubgraphMeanPathLengths():
    # Open Output File
    outFilename = 'tissue.ppi.mean_path_length'
    outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
    outFile = open(outFilepath, 'w')

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    # Iterate through all tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue
        graph = getTissueSubgraph(tissue)
        
        distances = getPathLengths(graph)
        distSum = 0
        numPaths = 0
        # Sum over all shortest paths
        for vertex in graph.vertices():
            distList = list(distances[vertex].a)
            for distance in distList:
                distSum = distSum + distance
                numPaths = numPaths + 1

        # This scheme double counts, so divide by two
        distSum = distSum / 2
        numPaths = numPaths / 2
        meanDist = distSum/float(numPaths)

        outFile.write('%s\t%f\n' % (tissue, meanDist))
        
# - - - - - DEGREE DISTRIBUTIONS - - - - - #

def getTissueSubgraphDegreeHistograms():

    AUGMENTED_TISSUE_LIST = ['global'] + FUNCTIONAL_TISSUE_LIST

    # Iterate through all tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        print tissue

        # Open Output File
        outFilename = '%s.ppi.degree_hist' % tissue
        outFilepath = PATH_TO_TISSUE_ANALYSIS + outFilename
        outFile = open(outFilepath, 'w')

        graph = getTissueSubgraph(tissue)

        buckets, counts = getDegreeHistogram(graph)

        numCounts = len(counts)

        for i in range(numCounts):
            outFile.write('%f\t%f\n' % (buckets[i], counts[i]))
        
