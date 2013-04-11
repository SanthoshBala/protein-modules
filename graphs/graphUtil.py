#! /usr/bin/python

# graphUtil.py
# Author: Santhosh Balasubramanian
# Created: February 16, 2013
# Last Modified: March 25, 2013


# Library Imports
from graph_tool.all import *

# Global Imports
from settings import *

# Utility Imports
from common.statistics import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


OUTLIER_THRESHOLD = 0.75


# - - - - - - - - - - GRAPH REPRESENTATION - - - - - - - - - - #


# getGeneSet: Returns set of all gene IDs representing vertices for <graph>.
def getGeneSet(graph):
    geneSet = set()
    geneIDProp = graph.vertex_properties['gene_id']

    # Iterate through Vertices
    for vertex in graph.vertices():
        geneSet.add( geneIDProp[vertex] )

    return geneSet

# getInteractionSet: Returns set of all interactions representing edges
# in <graph>. Represents interactions as a list of tuples between interactors.
def getInteractionSet(graph):
    interactionSet = set()
    geneIDProp = graph.vertex_properties['gene_id']

    # By default, the tuple (1, 2) != (2, 1). Moreover, you cannot have a
    # set of sets. Therefore, the only way to have a query-able edge set
    # is to have a set of tuples in which the tuples have a canonical 
    # ordering. The edge K<-->B is therefore represented as (B, K) (sorted).
    for edge in graph.edges():
        # Get Gene IDs for Interactors
        vSource = edge.source()
        vTarget = edge.target()
        vSourceID = geneIDProp[vSource]
        vTargetID = geneIDProp[vTarget]

        # Construct Sorted Tuple
        geneList = [vSourceID, vTargetID]
        geneList.sort()
        geneTuple = tuple(geneList)

        interactionSet.add(geneTuple)

    return interactionSet

# getIDMapForGraph: Returns dictionary for gene_id:vertex map for <graph>.
def getIDMapForGraph(graph):
    idVertexMap = {}

    # Get graph_tool's Internal gene_id:vertex Map
    geneIDProp = graph.vertex_properties['gene_id']

    # Iterate through Vertices
    for vertex in graph.vertices():
        geneID = geneIDProp[vertex]
        idVertexMap.update( { geneID : vertex } )
    
    return idVertexMap

# getAdjacencyMatrix: Returns adjacency matrix representing <graph>.
def getAdjacencyMatrix(graph):
    adjacencyMatrix = graph_tool.spectral.adjacency(graph)

    return adjacencyMatrix

# getLaplacianMatrix: Returns LaPlacian matrix representing <graph>.
def getLaplacianMatrix(graph):
    laplacianMatrix = graph_tool.spectral.laplacian(graph)

    return laplacianMatrix


# - - - - - - - - - - GRAPH SEARCH - - - - - - - - - - #


# getVertexListByProperty: Returns list of vertices in <graph> with 
# <vertexProperty> matching <propertyValue>.
def getVertexListByProperty(graph, vertexProperty, propertyValue):
    vertexList = graph_tool.util.find_vertex(graph, vertexProperty, 
                                             propertyValue)
    return vertexList

# getVertexByProperty: Returns 1 vertex in <graph> with <vertexProperty> 
# matching <propertyValue>. Returns None if no such vertex exists, and does
# not attempt to distinguish between multiple matching vertices. 
def getVertexByProperty(graph, vertexProperty, propertyValue):
    vertexList = getVertexListByProperty(graph, vertexProperty, propertyValue)

    if len(vertexList) == 0:
        return None
    else:
        return vertexList[0]

# getOrCreateVertex: Returns 1 vertex in <graph> with <vertexProperty>
# matching <propertyValue>. Creates one if no such vertex exists.
def getOrCreateVertex(graph, vertexProperty, propertyValue):
    # Try to Get Vertex
    vertex = getVertexByProperty(graph, vertexProperty, propertyValue)
    
    # Create Vertex
    if not vertex:
        vertex = graph.add_vertex()
        vertexProperty[vertex] = propertyValue

    return vertex

# getEdgeListByProperty: Returns list of edges in <graph> with <edgeProperty>
# matching <propertyValue>.
def getEdgeListByProperty(graph, edgeProperty, propertyValue):
    edgeList = graph_tool.util.find_edge(graph, edgeProperty, propertyValue)

    return edgeList

# getEdgeByProperty: Returns 1 edge in <graph> with <edgeProperty> matching
# <propertyValue>. Returns None if no such edge exists, and does not attempt
# to distinguish between multiple matching edges. Always returns 1st match.
def getEdgeByProperty(graph, edgeProperty, propertyValue):
    edgeList = getEdgeListByProperty(graph, edgeProperty, propertyValue)

    if len(edgeList) == 0:
        return None
    else:
        return edgeList[0]


# - - - - - - - - - - GLOBAL GRAPH PROPERTIES - - - - - - - - - - #


# getNumVertices: Returns # of vertices in <graph>.
def getNumVertices(graph):
    return graph.num_vertices()

# getNumEdges: Returns # of edges in <graph>.
def getNumEdges(graph):
    return graph.num_edges()

# getDensity: Returns density of <graph> -- E/[ V*(V-1)/2 ]
def getDensity(graph):
    numVertices = getNumVertices(graph)
    numerator = getNumEdges(graph)
    denominator = (numVertices)*(numVertices-1)/(2.0)
    density = numerator/denominator

    return density

# getDiameter: Returns diameter of <graph>.
def getDiameter(graph):
    # Compute Distance Histogram
    distanceBuckets, distanceCount = getShortestDistanceHistogram(graph)
    numBuckets = len(distanceCount)
    
    # Diameter = Value of Maximum Bucket
    iterList = range(numBuckets)
    iterList.reverse()

    # Iterate from End of List
    for i in iterList:
        # If Bucket is Non-Empty, Return Distance
        if distanceCount[i] > 0:
            return distanceBuckets[i]

# getMeanDegree: Returns average degree amongst vertices in <graph>.
def getMeanDegree(graph):
    meanDegree = graph_tool.stats.vertex_average(graph, 'out')

    return meanDegree

# getGlobalClusteringCoefficient: Computes average clustering coefficient 
# amongst all vertices in <graph>.
def getGlobalClusteringCoefficient(graph):
    gcc = graph_tool.clustering.global_clustering(graph)

    return gcc

# getConnectedness: Returns True if <graph> is connected.
def getConnectedness(graph):
    componentProperty, hist = graph_tool.topology.label_components(graph)

    return len(hist) == 1

# getBipartiteness: Returns True if <graph> is bipartite.
def getBipartiteness(graph):
    bipartiteness = graph_tool.topology.is_bipartite(graph)

    return bipartiteness

# getCentralPointDominance: Computes extent to which <graph> is clustered
# around a single central vertex.
def getCentralPointDominance(graph):
    vBetweenness, eBetweenness = graph_tool.centrality.betweenness(graph)
    cpd = graph_tool.centrality.central_point_dominance(graph, vBetweenness)

    return cpd

# computeMinSpanningTree: Returns minimum spanning tree of <graph>.
def computeMinSpanningTree(graph):
    minSpanningTree = graph_tool.topology.min_spanning_tree(graph)

    return minSpanningTree

# getVertexPropertyAverage: Returns average of <vertexProperty> in <graph>.
def getVertexPropertyAverage(graph, vertexProperty):
    vPropAvg = graph_tool.stats.vertex_average(graph, vertexProperty)

    return vPropAvg

# getEdgePropertyAverage: Returns average of <edgeProperty> in <graph>.
def getEdgePropertyAverage(graph, edgeProperty):
    ePropAvg = graph_tool.stats.edge_average(graph, edgeProperty)
    return ePropAvg


# - - - - - - - - - - LOCAL GRAPH PROPERTIES - - - - - - - - - - #


# getConnectedComponents: Returns vertexProperty labeling each vertex by
# which connected component it belongs to. Also returns <hist>, which
# indicates size of each connected component (indexed by component ID).
def getConnectedComponents(graph):
    componentProp, hist = graph_tool.topology.label_components(graph)

    return componentProp, hist

# computeBetweenness: Computes vertex betweenness and edge betweenness
# and adds to <graph> as vertex and edge properties, respectively
def computeBetweenness(graph):
    # Compute Betweenness
    vBetweenness, eBetweenness = graph_tool.centrality.betweenness(graph)

    # Add to Graph Properties
    graph.edge_properties['edge_betweenness'] = eBetweenness
    graph.vertex_properties['vertex_betweenness'] = vBetweenness

    return

# computeLocalClusteringCoefficients: Computes individual clustering
# coefficients for each vertex in <graph> and adds as vertex property.
def computeLocalClusteringCoefficients(graph):
    localClusters = graph_tool.clustering.local_clustering(graph)
    graph.vertex_properties['local_clustering_coefficient'] = localClusters

    return

# getVertexPropertyOutliers: Returns list of vertices in top X% in <graph>.
def getVertexPropertyOutliers(graph, vertexProperty):
    # Get Upper Bound
    valueArray = vertexProperty.get_array()
    valueList = valueArray.tolist()
    upperBound = max(valueList)

    # Get Lower Bound
    lowerBound = OUTLIER_THRESHOLD*upperBound

    # Get Outliers
    outliers = graph_tool.util.find_vertex_range(graph, vertexProperty,
                                                 (lowerBound, upperBound))

    return outliers

# getEdgePropertyOutliers: Returns list of edges in top X% in <graph>
def getEdgePropertyOutliers(graph, edgeProperty):
    # Get Upper Bound
    valueArray = edgeProperty.get_array()
    valueList = valueArray.tolist()
    upperBound = max(valueList)

    # Get Lower Bound
    lowerBound = OUTLIER_THRESHOLD*upperBound

    # Get Outliers
    outliers = graph_tool.util.find_edge_range(graph, edgeProperty,
                                               (lowerBound, upperBound))

    return outliers


# - - - - - - - - - - PROPERTY HISTOGRAMS - - - - - - - - - - #


# getVertexPropertyHistogram: Returns list of buckets and counts for those
# buckets, corresponding to <vertexProperty> for <graph>.
def getVertexPropertyHistogram(graph, vertexProperty):
    # Get Histogram
    histResult = graph_tool.stats.vertex_hist(graph, vertexProperty)

    # Convert Histogram to Lists
    countArray = histResult[0]
    bucketArray = histResult[1]
    countList = countArray.tolist()
    bucketList = bucketArray.tolist()

    return bucketList, countList

# getEdgePropertyHistogram: Returns list of buckets and counts for those
# buckets, corresponding to <edgeProperty> for <graph>.
def getEdgePropertyHistogram(graph, edgeProperty):
    # Get Histogram
    histResult = graph_tool.stats.edge_hist(graph, edgeProperty)

    # Convert Histogram to Lists
    countArray = histResult[0]
    bucketArray = histResult[1]
    countList = countArray.tolist()
    bucketList = countArray.tolist()

    return bucketList, countList

# getDegreeHistogram: Returns histogram of # vertices by degree.
def getDegreeHistogram(graph):
    # Get Degree Property Map
    degreeProp = graph.degree_property_map('out')
    bucketList, countList = getVertexPropertyHistogram(graph, degreeProp)

    return bucketList, countList

# getShortestDistanceHistogram: Returns histogram of path lengths.
def getShortestDistanceHistogram(graph):
    # Get Histogram
    histResult = graph_tool.stats.distance_histogram(graph)

    # Convert Histogram to Lists
    countArray = histResult[0]
    bucketArray = histResult[1]
    countList = countArray.tolist()
    bucketList = bucketArray.tolist()

    return bucketList, countList

# getPathLengthMatrix: Returns matrix of all pairwise path lengths in <graph>.
def getPathLengthMatrix(graph):
    pathLengths = graph_tool.topology.shortest_distance(graph)

    return pathLengths


# - - - - - - - - - - PAIRWISE GRAPH PROPERTIES - - - - - - - - - - #


# getIsomorphism: Returns True if <graphOne> and <graphTwo> are isomorphic.
def getIsomorphism(graphOne, graphTwo):
    isomorphism = graph_tool.topology.isomorphism(graphOne, graphTwo)
    return isomorphism

# getSimilarity: Returns True if <graphOne> and <graphTwo> are similar.
def getSimilarity(graphOne, graphTwo):
    similarity = graph_tool.topology.similarity(graphOne, graphTwo)
    return similarity

# getGeneSetUnion: Returns union of gene sets for <graphOne> and <graphTwo>.
def getGeneSetUnion(graphOne, graphTwo):
    # Get Gene Sets
    g1GeneSet = getGeneSet(graphOne)
    g2GeneSet = getGeneSet(graphTwo)

    # Get Union
    union = g1GeneSet.union(g2GeneSet)

    return union

# getGeneSetIntersection: Returns intersection of <graphOne> and <graphTwo> 
# gene sets.
def getGeneSetIntersection(graphOne, graphTwo):
    # Get Gene Sets
    g1GeneSet = getGeneSet(graphOne)
    g2GeneSet = getGeneSet(graphTwo)
    
    # Get Intersection
    intersection = g1GeneSet.intersection(g2GeneSet)

    return intersection

# getInteractionSetUnion: Returns union of interaction sets for <graphOne> and
# <graphTwo>.
def getInteractionSetUnion(graphOne, graphTwo):
    # Get Interaction Sets
    g1InteractionSet = getInteractionSet(graphOne)
    g2InteractionSet = getInteractionSet(graphTwo)
    
    # Get Union
    union = g1InteractionSet.union(g2InteractionSet)

    return union

# getInteractionSetIntersection: Returns intersection of interaction sets for
# <graphOne> and <graphTwo>.
def getInteractionSetIntersection(graphOne, graphTwo):
    # Get Interaction Sets
    g1InteractionSet = getInteractionSet(graphOne)
    g2InteractionSet = getInteractionSet(graphTwo)
    
    # Get Intersection
    intersection = g1InteractionSet.intersection(g2InteractionSet)

    return intersection

# getGraphUnion: Returns union of <graphOne> and <graphTwo. 
def getGraphUnion(graphOne, graphTwo):
    # Correlate with Gene ID
    g1IDProp = graphOne.vertex_properties['gene_id']
    g2IDProp = graphTwo.vertex_properties['gene_id']

    # Generate Union
    union = graph_tool.generation.graph_union(graphOne, graphTwo, 
                                              props = [(g1IDProp, g2IDProp)])

    return union

# getGraphIntersection: Returns intersection of <graphOne> and <graphTwo>.
def getGraphIntersection(graphOne, graphTwo):
    # Initialize Empty Intersection
    intersection = Graph(directed = False)
    intIDProp = intersection.new_vertex_property('string')
    intConfidenceProp = intersection.new_edge_property('double')
    intersection.vertex_propertes['gene_id'] = intIDProp
    intersection.edge_properties['confidence'] = intConfidenceProp
    

    # Get Interaction Intersection
    edgeIntersection = getInteractionSetIntersection(graphOne, graphTwo)

    # Track Vertices in Intersection
    intIDVertexMap = {}

    # Iterate through <edgeIntersection>
    for edge in edgeIntersection:
        sourceID = edge[0]
        targetID = edge[1]
        
        # Get <graphOne> Vertices
        g1IDProp = graphOne.vertex_properties['gene_id']
        g1Source = getVertexByProperty(graphOne, g1IDProp, sourceID)
        g1Target = getVertexByProperty(graphOne, g1IDProp, targetID)

        # Get <graphTwo> Vertices
        g2IDProp = graphTwo.vertex_properties['gene_id']
        g2Source = getVertexByProperty(graphTwo, g2IDProp, sourceID)
        g2Target = getVertexByProperty(graphTwo, g2IDProp, targetID)

        # Get or Create Vertices
        vSource = getOrCreateVertex(intersection, intIDProp, sourceID)
        vTarget = getOrCreateVertex(intersection, intIDProp, targetID)

        # Create Edge
        intEdge = intersection.add(vSource, vTarget)

        # Set Confidence Value
        g1ConfProp = graphOne.edge_properties['confidence']
        g2ConfProp = graphTwo.edge_properties['confidence']
        g1Edge = graphOne.edge(g1Source, g1Target)
        g2Edge = graphTwo.edge(g2Source, g2Target)
        g1EdgeConf = g1ConfProp[g1Edge]
        g2EdgeConf = g2ConfProp[g2Edge]
        meanConf = mean([g1EdgeConf, g2EdgeConf])
        intConfidenceProp[intEdge] = meanConf

    return intersection
