#! /usr/bin/python

# moduleTopology.py
# Author: Santhosh Balasubramanian
# Created: February 16, 2013
# Last Modified: March 27, 2013


# Python Imports
import os
import math
import copy
from multiprocessing import Process

# Library Imports
from graph_tool.all import *

# Global Imports
from settings import *

# Module Imports
from modules.moduleUtil import *
from modules.moduleDatabase import *

# Graph Imports
from graphs.graphIO import *
from graphs.graphUtil import *

# Ontology Imports
from ontology.geneOntology import *


# - - - - - - - - - - ALGORITHM - - - - - - - - - - #


# A "module topology" M is a tree used to describe the topology of
# an original graph G. To construct M, G is recursively clustered. Consider
# a single clustering of G, resulting in a list of modules, each of which
# contains a set of nodes in G. Each of these modules is represented as a 
# node in M, with children corresponding to the clustered nodes. 
# G is then converted to G' by collapsing all nodes belonging to 
# module m into a single node in G'. Recursive iteration results in a tree.
#
# There are two settings for the algorithm. INIT/NO_INIT describes whether
# or not M is initialized with all nodes in G. At the end of the algorithm,
# the root of M becomes the parent of any unclustered nodes. PROP/NO_PROP
# describes whether or not unclustered nodes are propagated to the next
# level of tree M. For example, if X is not clustered in round i - 1, then
# X is given a parent in M called X@i. Based on the settings of these 
# variables, three possible types of module topologies are possible, each
# with different properties.


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


PROPAGATE_GENES = False
INITIALIZE_MODULE_TOPOLOGY = False

MODULE_TOPOLOGY_HEIGHT_BASE_FILENAME = MODULE_TOPOLOGY_BASE_FILENAME + '.heights'
PATH_TO_SUPER_MODULES = PATH_TO_ANALYSIS + 'super_modules/'


# - - - - - - - - - - MODULE TOPOLOGY CREATION - - - - - - - - - - #


# createAllModuleTopologies: Creates module topologies for tissue
# ppi subgraphs in range.
def createAllModuleTopologies(minimum = 0, maximum = 85, 
                              initialize = INITIALIZE_MODULE_TOPOLOGY,
                              propagate = PROPAGATE_GENES):

    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST[minimum:maximum]:
        if PRINT_PROGRESS:
            print tissue

        # Open Tissue Network File
        interaction = DESIRED_INTERACTION_TYPE
        inFileName = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, interaction)
        inFilePath = PATH_TO_TISSUE_SUBGRAPHS + inFileName

        # Get Graph
        graph = readGraph(inFilePath)

        # Get Module Topology
        moduleTopology = getModuleTopology(graph)

        # Get Output Path
        if initialize:
            if propagate:
                outDir = PATH_TO_MODULE_TOPOLOGIES + 'init_prop/'
            else:
                outDir = PATH_TO_MODULE_TOPOLOGIES + 'init_no_prop/'
        else:
            outDir = PATH_TO_MODULE_TOPOLOGIES + 'no_init/'

        # Open Module Topology Output File
        outFileName = MODULE_TOPOLOGY_BASE_FILENAME % tissue
        outFilePath = outDir + outFileName

        # Write to File
        writeGraph(moduleTopology, outFilePath, canonical = True)
        
        # Write Heights to File
        heightProp = moduleTopology.vertex_properties['height']
        geneIDProp = moduleTopology.vertex_properties['gene_id']

        # Open Heights File
        heightsFilePath = outFilePath + '.heights'
        heightsFile = open(heightsFilePath, 'w')

        # Write Heights
        for vertex in moduleTopology.vertices():
            geneID = geneIDProp[vertex]
            height = heightProp[vertex]
            heightsFile.write('%s\t%d\n' % (geneID, height))

    return

# getModuleTopology: Returns module topology for <graph>.
def getModuleTopology(graph, initialize = INITIALIZE_MODULE_TOPOLOGY,
                      propagate = PROPAGATE_GENES):
    modTopology = Graph(directed=False)

    # Create Graph Properties
    modTopology.vertex_properties['height'] = modTopology.new_vertex_property('int')
    modTopology.vertex_properties['gene_id'] = modTopology.new_vertex_property('string')
    modTopology.vertex_properties['parent'] = modTopology.new_vertex_property('object')
    modTopology.edge_properties['confidence'] = modTopology.new_edge_property('double')

    modGeneIDProp = modTopology.vertex_properties['gene_id']
    modHeightProp = modTopology.vertex_properties['height']
    parentProp = modTopology.vertex_properties['parent']
    frontier = dict()
    
    graphGeneIDProp = graph.vertex_properties['gene_id']

    # Initialize Module Topology
    currentHeight = 0
    if initialize:
        for vertex in graph.vertices():
            newVertex = modTopology.add_vertex()
            if propagate:
                modGeneIDProp[newVertex] = '%s@%d' % (graphGeneIDProp[vertex], currentHeight)
            else:
                modGeneIDProp[newVertex] = '%s' % graphGeneIDProp[vertex]
            modHeightProp[newVertex] = currentHeight
            frontier.update( { graphGeneIDProp[vertex] : newVertex } )

    # Call Recursive moduleTopologyHelper
    moduleTopology = moduleTopologyHelper(graph, modTopology, frontier, currentHeight,
                                          initialize = initialize,
                                          propagate = propagate)
    
    return moduleTopology

# moduleTopologyHelper: <graph> is the current working version of the graph
# whose <moduleTopology> is being constructed. Thus, in lower levels of
# the recursion, <graph> will no longer represent the original graph.
# At each level, <frontier> contains a series of pointers to nodes in
# <moduleTopology> thaht were added in previous level of recursive call stack.
def moduleTopologyHelper(graph, moduleTopology, frontier, height,
                         initialize = INITIALIZE_MODULE_TOPOLOGY,
                         propagate = PROPAGATE_GENES):
    moduleTempFilename = PATH_TO_MODULES + '%s.modules.temp' % str(os.getpid())

    # Recursive Base Case 1 - Graph = 1 Node
    if graph.num_vertices() == 1:
        return moduleTopology

    # Run SPICi
    createGraphRawModuleFile(graph, moduleTempFilename)

    # Open Module Temp File
    moduleTempFile = open(moduleTempFilename, 'r')
    
    # Track <numModules> to Identify Base Case 2
    numModules = 0
    newFrontier = dict()
    parentProp = moduleTopology.vertex_properties['parent']
    geneIDProp = moduleTopology.vertex_properties['gene_id']
    heightProp = moduleTopology.vertex_properties['height']
    confidenceProp = moduleTopology.edge_properties['confidence']

    # Iterate through File
    for line in moduleTempFile:
        # Increment numModules
        numModules = numModules + 1

        # Construct Module Gene List
        moduleGeneNames = list()
        lineTabSplit = parseTabSeparatedLine(line)
        for item in lineTabSplit:
            if propagate:
                geneName = item.split('@')[0]
            else:
                geneName = item
            moduleGeneNames.append(geneName)

        # Sort moduleGeneNames
        moduleGeneNames.sort()

        # Collapse Nodes in <graph>
        graph = collapseModule(graph, moduleGeneNames)

        # Update module topology
        moduleTopology = updateModuleTopology(moduleTopology, moduleGeneNames,
                                              frontier, height,
                                              initialize=initialize,
                                              propagate=propagate)
        
    # Remove moduleTempFile
    moduleTempFile.close()
    os.remove(moduleTempFilename)

    # Recursive Base Case 2 - Graph has 0 Modules
    if numModules == 0:
        # Get Unclustered Nodes
        if not propagate:
            frontier = dict()
            for vertex in moduleTopology.vertices():
                if parentProp[vertex] == None:
                    geneName = geneIDProp[vertex]
                    frontier.update( { geneName : vertex } )

        moduleGeneNames = frontier.keys()
        moduleTopology = updateModuleTopology(moduleTopology, moduleGeneNames,
                                              frontier, height, 
                                              initialize=initialize,
                                              propagate=propagate)
        return moduleTopology
    else:
        # Set Parent Links for Unclustered Nodes
        if propagate:
            for geneName, geneNode in frontier.iteritems():
                if parentProp[geneNode] == None:
                    newNode = moduleTopology.add_vertex()
                    geneIDProp[newNode] = '%s@%d' % (geneName, height + 1)
                    heightProp[newNode] = height + 1
                    parentProp[geneNode] = newNode
                    newEdge = moduleTopology.add_edge(geneNode, newNode)
                    confidenceProp[newEdge] = 1.0
                else:
                    continue

    # Reset Frontier
    frontier = dict()
    if propagate:
        newFrontierNodes = getVertexListByProperty(moduleTopology, heightProp, 
                                                   height + 1)
        for node in newFrontierNodes:
            geneName = geneIDProp[node].split('@')[0]
            frontier.update( { geneName : node })
    else:
        if initialize:
            for vertex in moduleTopology.vertices():
                geneName = geneIDProp[vertex]
                frontier.update( { geneName : vertex } )
        else:
            newFrontierNodes = getVertexListByProperty(moduleTopology, 
                                                      heightProp, 
                                                      height + 1)
            for node in newFrontierNodes:
                geneName = geneIDProp[node]
                frontier.update( { geneName : node })

    # Return
    return moduleTopologyHelper(graph, moduleTopology, frontier, height + 1,
                                initialize=initialize, propagate=propagate)

# updateModuleTopology: Updates <moduleTopology> with the modules in 
# <moduleGeneNames>, which were derived from another graph.
def updateModuleTopology(moduleTopology, moduleGeneNames, frontier, height,
                         initialize = INITIALIZE_MODULE_TOPOLOGY,
                         propagate = PROPAGATE_GENES):
    # Get Graph Properties
    heightProp = moduleTopology.vertex_properties['height']
    geneIDProp = moduleTopology.vertex_properties['gene_id']
    parentProp = moduleTopology.vertex_properties['parent']
    confProp = moduleTopology.edge_properties['confidence']

    # Add New Module Node
    moduleNode = moduleTopology.add_vertex()
    heightProp[moduleNode] = height + 1
    if propagate:
        # Get or Create a Module ID
        moduleID = getModuleID(moduleGeneNames)
        if not moduleID:
            moduleID = 'mod:%s' % moduleGeneNames[0]
        geneIDProp[moduleNode] = '%s@%d' % (moduleID, height + 1)
    else:
        # Get or Create a Module ID
        moduleID = getModuleID(moduleGeneNames)
        if not moduleID:
            moduleID = 'mod:%s' % moduleGeneNames[0]
        geneIDProp[moduleNode] = moduleID

    # Add Edges Between all Genes in Module and Module Node
    for geneName in moduleGeneNames:
        if initialize:
            geneNode = frontier.get(geneName)
        else:
            geneNode = frontier.get(geneName)
            if not geneNode:
                geneNode = moduleTopology.add_vertex()
                geneIDProp[geneNode] = geneName
                heightProp[geneNode] = height
        newEdge = moduleTopology.add_edge(geneNode, moduleNode)
        confProp[newEdge] = 1.0
        parentProp[geneNode] = moduleNode

    return moduleTopology

# collapseModule: Collapses <moduleGeneNames> in <graph> to single node.
def collapseModule(graph, moduleGeneNames):
    # Get Gene ID Property
    geneIDProp = graph.vertex_properties['gene_id']
    edgeWeightProp = graph.edge_properties['confidence']

    moduleNode = graph.add_vertex()
    
    # Get or Create a Module ID
    moduleID = getModuleID(moduleGeneNames)
    if not moduleID:
        moduleID = 'mod:%s' % moduleGeneNames[0]

    geneIDProp[moduleNode] = moduleID

    # Iterate through Gene Names in module
    for geneName in moduleGeneNames:
        geneNode = getVertexByProperty(graph, geneIDProp, geneName)
        
        numNeighbors = 0
        edgeWeightSum = 0

        # Iterate over Edges to Compute Avg Edge Weight
        for edge in geneNode.all_edges():
            # Get the "other" Node
            if edge.source() == geneNode:
                otherNode = edge.target()
            else:
                otherNode = edge.source()

            # Verify that "otherNode" not in Module
            if geneIDProp[otherNode] in moduleGeneNames:
                continue

            # Else, add Edge Weight
            edgeWeightSum = edgeWeightSum + edgeWeightProp[edge]
            # Increment <numNeighbors>
            numNeighbors = numNeighbors + 1

        if numNeighbors == 0:
            continue

        # Get New Edge Weight
        avgEdgeWeight = edgeWeightSum/float(numNeighbors)

        # Iterate through <geneNode> Neighbors
        for neighbour in geneNode.all_neighbours():
            # If neighbour not in module, add moduleNode<->neighbour
            # However, need to determine the weight for edge.
            if geneIDProp[neighbour] not in moduleGeneNames:
                newEdge = graph.add_edge(neighbour, moduleNode)
                edgeWeightProp[newEdge] = avgEdgeWeight

            # Else, continue
            else:
                continue

    for geneName in moduleGeneNames:
        # Remove Module Component from Graph
        geneNode = getVertexByProperty(graph, geneIDProp, geneName)
        graph.remove_vertex(geneNode)
            
    return graph


# - - - - - - - - - - TOPOLOGY HEIGHT - - - - - - - - - - #        


# getTopologyHeight: Returns height of <topology>.
def getTopologyHeight(topology):
    heightProp = topology.vertex_properties['height']
    buckets, counts = getVertexPropertyHistogram(topology, heightProp)
    height = max(buckets)

    return height

# getAllTopologyHeights: Gets heights of all module topology trees.
def getAllTopologyHeights(initialize = INITIALIZE_MODULE_TOPOLOGY,
                          propagate = PROPAGATE_GENES):
    # Figure out Filename and Path
    if initialize:
        if propagate:
            inFilepath = PATH_TO_MODULE_TOPOLOGIES + 'init_prop/'
            outFilename = 'init.prop.topology.heights'
        else:
            inFilepath = PATH_TO_MODULE_TOPOLOGIES + 'init_no_prop/'
            outFilename = 'init.no_prop.topology.heights'
    else:
        inFilepath = PATH_TO_MODULE_TOPOLOGIES + 'no_init/'
        outFilename = 'no_init.topology.heights'

    outFile = open(PATH_TO_MODULE_TOPOLOGIES + outFilename, 'w')

    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue

        # Open Heights File
        inFileName = MODULE_TOPOLOGY_HEIGHT_BASE_FILENAME % tissue
        inFile = open(inFileName, 'r')

        # Get Height
        height = 0
        for line in inFile:
            lineFields = parseTabSeparatedLine(line)
            if int(lineFields[1]) > height:
                height = int(lineFields[1])

        outFile.write('%s\t%d\n' % (tissue, height))

    # Close Files
    outFile.close()

    return
    

# - - - - - - - - - - TOPOLOGY BRANCHING SERIES - - - - - - - - - - #


# getTopologyBranchingSeries: For each height, returns average branching #.
def getTopologyBranchingSeries(topology):
    height = getTopologyheight(topology)
    heightProp = topology.vertex_properties['height']

    branchingSeries = list()

    # Iterate through Tree Levels
    for i in range(height):
        vertexSet = getVertexListByProperty(topology, heightProp, i)
        numVertices = len(vertexSet)

        # Iterate through Level Vertices
        degreeSum = 0
        for vertex in vertexSet:
            degreeSum = degreeSum + vertex.out_degree()
        
        # Subtract Parents from Degree Computation
        branchingSum = degreeSum - numVertices
        avgBranchingNum = branchingSum/float(numVertices)
        
        branchingSeries.append(avgBranchingNum)
        
    return branchingSeries

# getAllTopologyBranchingSeries: Gets branching series for all topologies.
def getAllTopologyBranchingSeries(initialize = INITIALIZE_MODULE_TOPOLOGY,
                                  propagate = PROPAGATE_GENES):

    # Figure out Filename and Path
    if initialize:
        if propagate:
            inFilepath = PATH_TO_MODULE_TOPOLOGIES + 'init_prop/'
            outFilename = 'init.prop.topology.branching'
        else:
            inFilepath = PATH_TO_MODULE_TOPOLOGIES + 'init_no_prop/'
            outFilename = 'init.no_prop.topology.branching'
    else:
        inFilepath = PATH_TO_MODULE_TOPOLOGIES + 'no_init/'
        outFilename = 'no_init.topology.branching'

    outFile = open(PATH_TO_MODULE_TOPOLOGIES + outFilename, 'w')

    # Iterate through Tissues
    for tissue in AUGMENTED_TISSUE_LIST:
        if PRINT_PROGRESS:
            print tissue

        graph = getTissueSubgraph()
        topology = getModuleTopology(graph)
        
        series = getTopologyBranchingSeries(topology)
        outFile.write('%s\t%s\n' % (tissue, str(series)))

    # Close File
    outFile.close()

    return


# - - - - - - - - - - MODULE TOPOLOGY ANALYSIS - - - - - - - - - - #


# createOntologyExpansionFiles
def createOntologyExpansionFiles():
    
    # Iterate through All Tissues
    for tissue in FUNCTIONAL_TISSUE_LIST:
        print tissue
        # Open Output File
        functionFilename = FUNCTION_EXPANSION_BASE_FILENAME % tissue
        processFilename = PROCESS_EXPANSION_BASE_FILENAME % tissue

        functionFile = open(PATH_TO_FUNCTION_EXPANSION + functionFilename,'w')
        processFile = open(PATH_TO_PROCESS_EXPANSION + processFilename, 'w')

        # Get Tissue Subgraph
        graph = getTissueSubgraph(tissue)

        # Get Module Topology
        topology = getModuleTopology(graph)
        geneIDProp = topology.vertex_properties['gene_id']

        # Get Function and Process Annotations
        functions, processes = getOntologyAnnotation(topology)

        # Get 0 Height Nodes
        heightProp = topology.vertex_properties['height']
        zeroHeightNodes = getVertexListByProperty(topology, heightProp, 0)

        # For each leaf, add line to output files
        parentProp = topology.vertex_properties['parent']
        for leaf in zeroHeightNodes:
            leafGeneID = geneIDProp[leaf]

            # Initialize Function and Process Series
            functionSeries = []
            processSeries = []
            leafFunctions = functions.get(leafGeneID)
            leafProcesses = processes.get(leafGeneID)
            if leafFunctions == None:
                continue
            if leafProcesses == None:
                continue
            functionSeries.append(len(leafFunctions))
            processSeries.append(len(leafProcesses))
            
            # Iterate Recursively to Root, filling in function/process series
            node = leaf
            while parentProp[node] != None:
                node = parentProp[node]
                nodeID = geneIDProp[node]
                functionSeries.append( len(functions.get(nodeID)) )
                processSeries.append( len(processes.get(nodeID)) )

            # Write results to files
            functionFile.write('%s\t%s\n' % (leafGeneID, str(functionSeries)))
            processFile.write('%s\t%s\n' % (leafGeneID, str(processSeries)))
    
    functionFile.close()
    processFile.close()

    return

# getTopologyAnnotation: Note that this function assumes the topology
# was formed by INITIALIZATION without PROPAGATION.
def getTopologyAnnotation(topology):
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    geneIDMap = db.geneIDMap

    # Get Gene ID Map
    geneIDProp = topology.vertex_properties['gene_id']

    # Get parent of a given gene
    parentProp = topology.vertex_properties['parent']

    # Initialize topologyAnnotation ( gene_id --> set() )
    functionAnnotation = {}
    processAnnotation = {}
    componentAnnotation = {}
    
    # Initialize topology Annotation with nodes at height 0
    height = getTopologyHeight(topology)
    heightProp = topology.vertex_properties['height']
    heightZeroNodes = getVertexListByProperty(topology, heightProp, 0)

    for node in heightZeroNodes:
        geneID = geneIDProp[node]
        
        # Convert 'ENSG' to Entrez if Necessary
        if 'ENSG' in geneID:
            funcs, procs, comps = getGeneAnnotation(geneID, 
                                                    nomenclature = 'ensembl', 
                                                    combined = False)
        else:
            funcs, procs, comps = getGeneAnnotation(geneID, combined = False)

        if not funcs or not procs or not comps:
            continue

        # Get Sets
        funcTuples = []
        procTuples = []
        compTuples = []
        for func in funcs:
            funcTuples.append(tuple(func))
        for proc in procs:
            procTuples.append(tuple(proc))
        for comp in comps:
            compTuples.append(tuple(comp))
        funcSet = set(funcTuples)
        procSet= set(procTuples)
        compSet = set(compTuples)

        functionAnnotation.update( { geneID : funcSet } )
        processAnnotation.update( { geneID : procSet } )
        componentAnnotation.update( { geneID : compSet } )

    # Iterate through height 1 to height
    for i in range(1, height):
        levelNodes = getVertexListByProperty(topology, heightProp, i)

        # Iterate through all nodes at this level
        for node in levelNodes:
            geneID = geneIDProp[node]
            
            nodeFunctions = set()
            nodeProcesses = set()
            nodeComponents = set()

            # Get children
            neighbors = node.all_neighbours()
            children = []
            parent = parentProp[node]
            for neighbor in neighbors:
                if neighbor == parent:
                    continue
                children.append(geneIDProp[neighbor])
                
            # Get union of annotations for children
            for child in children:
                childFunctions = functionAnnotation.get(child)
                childProcesses = processAnnotation.get(child)
                childComponents = componentAnnotation.get(child)

                # We don't have annotations for all nodes, so
                # handle this case gracefully
                if childFunctions == None:
                    continue
                if childProcesses == None:
                    continue
                if childComponents == None:
                    continue

                if len(nodeFunctions) == 0:
                    nodeFunctions = copy.copy(childFunctions)
                    nodeProcesses = copy.copy(childProcesses)
                    nodeComponents = copy.copy(childComponents)
                else:
                    nodeFunctions = nodeFunctions.intersection(childFunctions)
                    nodeProcesses = nodeProcesses.intersection(childProcesses)
                    nodeComponents = nodeComponents.intersection(childComponents)
            
            # Add annotations to dictionaries
            functionAnnotation.update( { geneID : nodeFunctions } )
            processAnnotation.update( { geneID : nodeProcesses } )
            componentAnnotation.update( { geneID : nodeComponents } )

    return functionAnnotation, processAnnotation, componentAnnotation

# getSuperModules
def getSuperModules(tissue):
    # Get Tissue Subgraph
    graph = getTissueSubgraph(tissue)

    # Get Module Topology
    topology = getModuleTopology(graph)
    heightProp = topology.vertex_properties['height']
    idProp = topology.vertex_properties['gene_id']

    # Get Annotations
    functions, processes, components = getTopologyAnnotation(topology)

    # Get Height
    height = getTopologyHeight(topology)

    # Open Output File
    outFilePath = PATH_TO_SUPER_MODULES
    outFileName = '%s.super.modules' % tissue
    outFile = open(outFilePath + outFileName, 'w')

    # Iterate through Possible Super Modules
    for i in range(2, height):
        vertexList = getVertexListByProperty(topology, heightProp, i)

        # Iterate through Vertices at this Height
        for vertex in vertexList:
            moduleID = idProp[vertex]
            outFile.write('%s\n' % moduleID)
            vertexFuncs = functions.get(moduleID)
            for vFunc in vertexFuncs:
                outFile.write('\t%s\n' % str(vFunc))
            vertexProcs = processes.get(moduleID)
            for vProc in vertexProcs:
                outFile.write('\t%s\n' % str(vProc))
            vertexComps = components.get(moduleID)
            for vComp in vertexComps:
                outFile.write('\t%s\n' % str(vComp))

    outFile.close()
    
    return

