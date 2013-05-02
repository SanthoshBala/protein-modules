#! /usr/bin/python

from settings import *
from common import *
from graph_tool.all import *
from graphUtil import *
from moduleUtil import *

GRAPH_PROPERTIES = [
    'numVertices',
    'numEdges',
    'density',
    'connectedness',
    'centralDominance',
    'globalClusteringCoefficient',
    'averageDegree',
    'diameter',
    'bipartiteness',
    ]

GRAPH_PAIRWISE_PROPERTIES = [
    'vertex',
    'edge',
    'module',    
]

GENE_VARIABILITY_MEASURES = [
    'degree',
    'neighborhood',
    'module_partnership',
    'module_interaction',
    'local_cluster_coeff',
    'vertex_betweenness'
]

# analyzeTissueSubgraphsPairwise
def analyzeTissueSubgraphsPairwise():
    numTissues = len(GERM_LAYER_TISSUE_LIST)

    # Initialize Union Matrix Dict
    unionMatrixDict = {}
    for prop in GRAPH_PAIRWISE_PROPERTIES:
        matrix = constructSquareMatrix(numTissues)
        unionMatrixDict.update( { prop : matrix } )

    # Initialize Intersection Matrix Dict
    intMatrixDict = {}
    for prop in GRAPH_PAIRWISE_PROPERTIES:
        matrix = constructSquareMatrix(numTissues)
        intMatrixDict.update( { prop : matrix } )

    # Initialize Jaccard Matrix Dict
    jaccardMatrixDict = {}
    for prop in GRAPH_PAIRWISE_PROPERTIES:
        matrix = constructSquareMatrix(numTissues)
        jaccardMatrixDict.update( { prop : matrix } )

    # Iterate through all tissues
    for i in range(numTissues):
        print i
        tissueA = GERM_LAYER_TISSUE_LIST[i]
        graphA = getTissueSubgraph(tissueA)
        
        # Get Vertex Set
        verticesA = getVertexSet(graphA)

        # Get Edge Set
        edgesA = getEdgeSet(graphA)

        # Get Module Set
        modulesA = getTissueModuleIDSet(tissueA)

        # Iterate through all tissue pairs
        for j in range(i, numTissues):
            print '\t%d' % j
            tissueB = GERM_LAYER_TISSUE_LIST[j]
            graphB = getTissueSubgraph(tissueB)
        
            # Get Vertex Set
            verticesB = getVertexSet(graphB)
            # Get Edge Set
            edgesB = getEdgeSet(graphB)
            # Get Module Set
            modulesB = getTissueModuleIDSet(tissueB)
            
            # Get Vertex Union
            vertexUnion = verticesA.union(verticesB)
            # Get Vertex Intersection
            vertexInt = verticesA.intersection(verticesB)
            # Get Edge Union
            edgeUnion = edgesA.union(edgesB)
            # Get Edge Intersection
            edgeInt = edgesA.intersection(edgesB)
            # Get Module Union
            moduleUnion = modulesA.union(modulesB)
            # Get Module Intersection
            moduleInt = modulesA.intersection(modulesB)

            # Fill in Vertex Matrices
            vertUnionMatrix = unionMatrixDict.get('vertex')
            unionSize = float(len(vertexUnion))
            vertUnionMatrix[i][j] = unionSize
            vertUnionMatrix[j][i] = unionSize

            vertIntMatrix = intMatrixDict.get('vertex')
            intSize = float(len(vertexInt))
            if i == j:
                vertIntMatrix[i][j] = 0
                vertIntMatrix[j][i] = 0
            else:
                vertIntMatrix[i][j] = intSize
                vertIntMatrix[j][i] = intSize

            vertJaccardMatrix = jaccardMatrixDict.get('vertex')
            jaccard = intSize/unionSize
            if i == j:
                vertJaccardMatrix[i][j] = 0.0
                vertJaccardMatrix[j][i] = 0.0
            else:
                vertJaccardMatrix[i][j] = jaccard
                vertJaccardMatrix[j][i] = jaccard
            
            # Fill in Edge Matrices
            edgeUnionMatrix = unionMatrixDict.get('edge')
            unionSize = float(len(edgeUnion))
            edgeUnionMatrix[i][j] = unionSize
            edgeUnionMatrix[j][i] = unionSize

            edgeIntMatrix = intMatrixDict.get('edge')
            intSize = float(len(edgeInt))
            if i == j:
                edgeIntMatrix[i][j] = 0
                edgeIntMatrix[j][i] = 0
            else:
                edgeIntMatrix[i][j] = intSize
                edgeIntMatrix[j][i] = intSize

            edgeJaccardMatrix = jaccardMatrixDict.get('edge')
            jaccard = intSize/unionSize
            if i == j:
                edgeJaccardMatrix[i][j] = 0.0
                edgeJaccardMatrix[j][i] = 0.0
            else:
                edgeJaccardMatrix[i][j] = jaccard
                edgeJaccardMatrix[j][i] = jaccard
            
            # Fill in Module Matrices
            modUnionMatrix = unionMatrixDict.get('module')
            unionSize = float(len(moduleUnion))
            modUnionMatrix[i][j] = unionSize
            modUnionMatrix[j][i] = unionSize

            modIntMatrix = intMatrixDict.get('module')
            intSize = float(len(moduleInt))
            if i == j:
                modIntMatrix[i][j] = 0
                modIntMatrix[j][i] = 0
            else:
                modIntMatrix[i][j] = intSize
                modIntMatrix[j][i] = intSize

            modJaccardMatrix = jaccardMatrixDict.get('module')
            jaccard = intSize/unionSize
            if i == j:
                modJaccardMatrix[i][j] = 0.0
                modJaccardMatrix[j][i] = 0.0
            else:
                modJaccardMatrix[i][j] = jaccard
                modJaccardMatrix[j][i] = jaccard
            
    # Write Matrices to File
    for prop in GRAPH_PAIRWISE_PROPERTIES:
        intFilename = '%s.germ.intersection.pairwise.analysis' % prop
        intersectionFile = open(PATH_TO_PAIRWISE_ANALYSIS + intFilename, 'w')
        intMatrix = intMatrixDict.get(prop)
        writeMatrixToFile(intMatrix, intersectionFile)
        intersectionFile.close()

        unionFilename = '%s.germ.union.pairwise.analysis' % prop
        unionFile = open(PATH_TO_PAIRWISE_ANALYSIS + unionFilename, 'w')
        unionMatrix = unionMatrixDict.get(prop)
        writeMatrixToFile(unionMatrix, unionFile)
        unionFile.close()

        jaccardFilename = '%s.germ.jaccard.pairwise.analysis' % prop
        jaccardFile = open(PATH_TO_PAIRWISE_ANALYSIS + jaccardFilename, 'w')
        jaccardMatrix = jaccardMatrixDict.get(prop)
        writeMatrixToFile(jaccardMatrix, jaccardFile)
        jaccardFile.close()

# analyzeTissueSubgraphIntersections
def analyzeTissueSubgraphIntersections():
    numTissues = len(CANONICAL_TISSUE_LIST)

    # Initialize Union Matrix Dict
    unionMatrixDict = {}
    for prop in GRAPH_PROPERTIES:
        matrix = constructSquareMatrix(numTissues)
        unionMatrixDict.update( { prop : matrix } )

    # Initialize Intersection Matrix Dict
    intMatrixDict = {}
    for prop in GRAPH_PROPERTIES:
        matrix = constructSquareMatrix(numTissues)
        intMatrixDict.update( { prop : matrix } )

    # Similarity
    simMatrix = numTissues*[None]
    for i in range(len(simMatrix)):
        simMatrix[i] = numTissues*[None]
        
    # Isomorphism
    isoMatrix = numTissues*[None]
    for i in range(len(isoMatrix)):
        isoMatrix[i] = numTissues*[None]

    # First Tissue
    for i in range(numTissues):
        tissueA = CANONICAL_TISSUE_LIST[i]
        # Second Tissue
        for j in range(i, numTissues):
            tissueB = CANONICAL_TISSUE_LIST[j]
            
            # Get Graphs
            tissueAGraph, reportA = analyzeGenemaniaTissueSubgraph(tissueA)
            tissueBGraph, reportB = analyzeGenemaniaTissueSubgraph(tissueB)

            # Get Union
            union = getGraphUnion(tissueAGraph, tissueBGraph)
            union, unionReport = analyzeGraphToolObject(union)
            for key, value in unionReport.iteritems():
                if type(value) == list:
                    continue
                if key == 'averageDegree':
                    value = value[0]
                if key == 'globalClusteringCoefficient':
                    value = value[0]
                keyMatrix = unionMatrixDict[key]
                keyMatrix[i][j] = value
                keyMatrix[j][i] = value

            # Get Intersection
            intersection = getGraphIntersection(tissueAGraph, tissueBGraph)
            intersection, intReport = analyzeGraphToolObject(intersection)
            for key, value in intReport.iteritems():
                if type(value) == list:
                    continue
                if key == 'averageDegree':
                    value = value[0]
                if key == 'globalClusteringCoefficient':
                    value = value[0]
                keyMatrix = intMatrixDict[key]
                keyMatrix[i][j] = value
                keyMatrix[j][i] = value

            # Get Similarity
            similarity = getSimilarity(tissueAGraph, tissueBGraph)
            simMatrix[i][j] = similarity
            simMatrix[j][i] = similarity
            
            # Get Isomorphism
            isomorphism = getIsomorphism(tissueAGraph, tissueBGraph)
            isoMatrix[i][j] = isomorphism
            isoMatrix[j][i] = isomorphism
        print tissueA
    
    # Write Matrices to Files
    unionFile = open(PATH_TO_TISSUE_ANALYSIS + 'union.pairwise.analysis', 'w')
    unionFile.write('Analyzing Tissue Pairwise Union...\n\n')
    for prop, propMatrix in unionMatrixDict.iteritems():
        unionFile.write('%s\n' % prop)
        for row in propMatrix:
            for column in row:
                unionFile.write('\t%s' % str(column))
            unionFile.write('\n')
        unionFile.write('\n')
    unionFile.close()
    
    intFile = open(PATH_TO_TISSUE_ANALYSIS + 'int.pairwise.analysis', 'w')
    intFile.write('Analyzing Tissue Pairwise Intersection...\n\n')
    for prop, propMatrix in intMatrixDict.iteritems():
        intFile.write('%s\n' % prop)
        for row in propMatrix:
            for column in row:
                intFile.write('\t%s' % str(column))
            intFile.write('\n')
        intFile.write('\n')
    intFile.close()
 
    simFile = open(PATH_TO_TISSUE_ANALYSIS + 'similarity.pairwise.analysis', 'w')
    simFile.write('Analyzing Tissue Pairwise Similarity...\n\n')
    for row in simMatrix:
        for column in row:
            simFile.write('\t%f' % column)
        simFile.write('\n')
    simFile.write('\n')
    simFile.close()

    isoFile = open(PATH_TO_TISSUE_ANALYSIS + 'isomorphism.pairwise.analysis', 'w')
    isoFile.write('Analyzing Tissue Pairwise Isomorphism...\n\n')
    for row in isoMatrix:
        isoFile.write('\t')
        for column in row:
            isoFile.write('%f\t' % column)
        isoFile.write('\n')
    isoFile.write('\t')
    isoFile.close()

    return

# generateGlobalIntersectionGraph
def generateGlobalIntersectionGraph():
    intersection = Graph()
    
    # Set intersection = first graph
    copyGraph = True

    # Iterate through all tissues
    for tissue in CANONICAL_TISSUE_LIST:
        inFilename = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, DESIRED_INTERACTION_TYPE)
        tissueGraph = readGenemaniaGraph(PATH_TO_TISSUE_SUBGRAPHS + inFilename)

        if copyGraph:
            intersection = Graph(tissueGraph)
            copyGraph = False
            print tissue
            continue

        intersection = getGraphIntersection(intersection, tissueGraph)
        print tissue

    # Write out graph to disk
    outFilename = 'global.intersection.network'
    outFilePath = PATH_TO_TISSUE_SUBGRAPHS + outFilename
    outFile = open(outFilePath, 'w')

    # Iterate through all edges in graph
    outFile.write('Gene_A\tGene_B\tWeight\n')
    geneIDProp = intersection.vertex_properties['gene_id']
    confProp = intersection.edge_properties['confidence']
    for edge in intersection.edges():
        geneANode = edge.source()
        geneAName = geneIDProp[geneANode]
        geneBNode = edge.target()
        geneBName = geneIDProp[geneBNode]
        edgeConf = confProp[edge]
        outFile.write('%s\t%s\t%f\n' % (geneAName, geneBName, edgeConf))
        
    outFile.close()
    return intersection

# analyzeGlobalIntersectionGraph
def analyzeGlobalIntersectionGraph():
    inFilename = 'global.intersection.network'
    inFilePath = PATH_TO_TISSUE_SUBGRAPHS + inFilename

    graph = readGenemaniaGraph(inFilePath)
    tissueGraph, tissueReport = analyzeGraphToolObject(graph)

    outFile = open(PATH_TO_TISSUE_ANALYSIS + 'global.intersection.analysis', 
                   'w')

    for key, value in tissueReport.iteritems():
        outFile.write('%s\t%s\n' % (key, str(value)))

    geneIDProp = tissueGraph.vertex_properties['gene_id']

    # Write Degree Outliers
    degreeOutliers = getDegreeOutliers(tissueGraph)
    degreeProp = tissueGraph.degree_property_map('out')
    outFile.write('Degree Outliers\n')
    for outlier in degreeOutliers:
        outFile.write('\t%s\t%f\n' % (geneIDProp[outlier], 
                                      degreeProp[outlier]))
            
    # Write Vertex Betweenness Outliers
    vBetweennessOutliers = getVertexBetweennessOutliers(tissueGraph)
    vBetweennessProp = tissueGraph.vertex_properties['vertex_betweenness']
    outFile.write('Vertex Betweenness Outliers\n')
    for outlier in vBetweennessOutliers:
        outFile.write('\t%s\t%f\n' % (geneIDProp[outlier], 
                                      vBetweennessProp[outlier]))

# - - - - - HOUSEKEEPING GENES DB - - - - - #

# constructHousekeepingGenes
def constructHousekeepingGenesDB():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    housekeepingDB = db.housekeepingGenes

    housekeepingDict = {}
    housekeepingSet = set()

    # Initialize Housekeeping Set
    initSet = True

    # Iterate through all tissues
    numTissues = len(FUNCTIONAL_TISSUE_LIST)
    for i in range(numTissues):
        tissue = FUNCTIONAL_TISSUE_LIST[i]
        tissueGraph = getTissueSubgraph(tissue)
        tissueGraphGenes = getVertexSet(tissueGraph)

        if initSet:
            housekeepingSet = tissueGraphGenes
            initSet = False
            print tissue
            continue

        housekeepingSet = housekeepingSet.intersection(tissueGraphGenes)
        print tissue

    # Check intersection
    for tissue in FUNCTIONAL_TISSUE_LIST:
        tissueGraph = getTissueSubgraph(tissue)
        geneSet = getVertexSet(tissueGraph)
        for gene in housekeepingSet:
            if gene not in geneSet:
                print tissue, gene

    for gene in housekeepingSet:
        housekeepingDB.insert( { 'gene_id' : gene,
                                 'degree' : [],
                                 'neighbors' : [],
                                 'vertex_betweenness' : [],
                                 'local_clustering_coefficient' : [],
                                 'module_participation' : [],
                                 'module_identity' : [],
                                 'module_interaction' : [],                                 
                                 } )
    cli.close()

def fillHousekeepingGenesDB():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    housekeepingDB = db.housekeepingGenes

    # Iterate through all tissues
    numTissues = len(FUNCTIONAL_TISSUE_LIST)
    for i in range(numTissues):
        # Get graph
        tissue = FUNCTIONAL_TISSUE_LIST[i]
        print tissue
        tissueGraph = getTissueSubgraph(tissue)
        geneIDMap = tissueGraph.vertex_properties['gene_id']

        graphIDMap = getIDMapForGraph(tissueGraph)

        computeLocalClusteringCoefficients(tissueGraph)
        lccMap = tissueGraph.vertex_properties['local_clustering_coefficient']
        
        computeVertexBetweenness(tissueGraph)
        vbMap = tissueGraph.vertex_properties['vertex_betweenness']

        
        modParticipation = getModularizedGeneSet(tissue)


        # modIdentity
#        modInteraction = 

        # Iterate through all genes
        for gene_id, vertex in graphIDMap.iteritems():
            record = housekeepingDB.find_one( { 'gene_id' : gene_id } )

            if not record:
                continue
 
            # Degree
            degree = vertex.out_degree()
            record['degree'].append(degree)

            # LCC
            lcc = lccMap[vertex]
            record['local_clustering_coefficient'].append(lcc)

            # VB
            vb = vbMap[vertex]
            record['vertex_betweenness'].append(vb)

            # Module Participation
            if gene_id in modParticipation:
                record['module_participation'].append(True)
            else:
                record['module_participation'].append(False)

            # Module Identity

            # Module Interactors

            

            # Neighbors
            neighbors = list()
            for neighbor in vertex.out_neighbours():
                neighbors.append(geneIDMap[neighbor])
            record['neighbors'].append(neighbors)

            housekeepingDB.save(record)


# analyzeHousekeepingGenesAcrossTissues
def analyzeHousekeepingGenesAcrossTissues():
    intFilename = 'global.intersection.network'
    intFilePath = PATH_TO_TISSUE_SUBGRAPHS + intFilename

    outFilePath = PATH_TO_TISSUE_ANALYSIS + 'housekeeping.analysis'
    outFile = open(outFilePath, 'w')
    outFile.write('Analyzing Housekeeping Genes...\n\n')

    intersectionGraph = readGenemaniaGraph(intFilePath)
    intGeneID = intersectionGraph.vertex_properties['gene_id']
    dataDict = {}
    for geneVertex in intersectionGraph.vertices():
        geneName = intGeneID[geneVertex]
        dataDict.update( { geneName : {
                    'degree' : [],
                    'vb' : [],
                    'lcc' : [],
                    }
                           }
                         )

    for tissue in FUNCTIONAL_TISSUE_LIST:
        inFilename = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, 'ppi')
        inFilePath = PATH_TO_TISSUE_SUBGRAPHS + inFilename
        tissueGraph = readGenemaniaGraph(inFilePath)
        computeVertexBetweenness(tissueGraph)
        computeLocalClusteringCoefficients(tissueGraph)
        tissueGeneID = tissueGraph.vertex_properties['gene_id']
        degreeProp = tissueGraph.degree_property_map('out')
        lccProp = tissueGraph.vertex_properties['local_clustering_coefficient']
        vbProp = tissueGraph.vertex_properties['vertex_betweenness']
        # Iterate through vertices of intersection graph
        for intGeneVertex in intersectionGraph.vertices():
            intGeneName = intGeneID[intGeneVertex]
            
            
            tissueVertex = getVertexByProperty(tissueGraph, tissueGeneID, 
                                               intGeneName)[0]
            dataDict[intGeneName]['degree'].append(degreeProp[tissueVertex])
            dataDict[intGeneName]['lcc'].append(lccProp[tissueVertex])
            dataDict[intGeneName]['vb'].append(vbProp[tissueVertex])
        print tissue

    for gene, data in dataDict.iteritems():
        outFile.write('%s\n' % gene)
        outFile.write('\tdegree\t%s\n' % str(data['degree']))
        outFile.write('\tvb\t%s\n' % str(data['vb']))
        outFile.write('\tlcc\t%s\n' % str(data['lcc']))

# analyzeTissueSubgraphsIndividually
def analyzeTissueSubgraphsIndividually():
    degreeOutlierDict = {}
    vbOutlierDict = {}

    # Iterate through all tissues
    for tissue in CANONICAL_TISSUE_LIST:
        tissueGraph, tissueReport = analyzeGenemaniaTissueSubgraph(tissue)
        outFilename = '%s.%s.subgraph.analysis' % (tissue, 
                                                   DESIRED_INTERACTION_TYPE)
        outFile = open(PATH_TO_TISSUE_ANALYSIS + outFilename, 'w')

        outFile.write('Analyzing %s %s Network...\n\n' % 
                      (tissue, DESIRED_INTERACTION_TYPE))
        
        for key, value in tissueReport.iteritems():
            outFile.write('%s\t%s\n' % (key, str(value)))

        geneIDProp = tissueGraph.vertex_properties['gene_id']
        
        # Write Degree Outliers
#        degreeOutliers = getDegreeOutliers(tissueGraph)
#        degreeProp = tissueGraph.degree_property_map('out')
        outFile.write('Degree Outliers\n')
#        for outlier in degreeOutliers:
#            outFile.write('\t%s\t%f\n' % (geneIDProp[outlier], 
#                                          degreeProp[outlier]))
#            if geneIDProp[outlier] in degreeOutlierDict:
#                degreeOutlierDict[geneIDProp[outlier]].append(tissue)
#            else:
#                degreeOutlierDict.update( { geneIDProp[outlier] : [tissue] } )

        # Write Vertex Betweenness Outliers
#        vBetweennessOutliers = getVertexBetweennessOutliers(tissueGraph)
#        vBetweennessProp = tissueGraph.vertex_properties['vertex_betweenness']
        outFile.write('Vertex Betweenness Outliers\n')
#        for outlier in vBetweennessOutliers:
#            outFile.write('\t%s\t%f\n' % (geneIDProp[outlier], 
#                                          vBetweennessProp[outlier]))
#            if geneIDProp[outlier] in vbOutlierDict:
#                vbOutlierDict[geneIDProp[outlier]].append(tissue)
#            else:
#                vbOutlierDict.update( { geneIDProp[outlier] : [tissue] } )

   
        # Write Edge Betweenness Outliers
#        eBetweennessOutliers = getEdgeBetweennessOutliers(tissueGraph)
#        eBetweennessProp = tissueGraph.edge_properties['edge_betweenness']
        outFile.write('Edge Betweenness Outliers\n')
#        for outlier in eBetweennessOutliers:
#            orig = outlier.source()
#            dest = outlier.target()
#            outFile.write('\t%s\t%s\t%f\n' % (geneIDProp[orig], 
#                                              geneIDProp[dest], 
#                                              eBetweennessProp[outlier]))

        # Write Local Clustering Coefficient Outliers
#        lccOutliers = getLCCOutliers(tissueGraph)
#        lccProp = tissueGraph.vertex_properties['local_clustering_coefficient']
        outFile.write('Local Clustering Coefficient Outliers\n')
#        for outlier in lccOutliers:
#            outFile.write('\t%s\t%f\n' % (geneIDProp[outlier], 
#                                          lccProp[outlier]))

        print tissue

    # Write outliers file
    f = open(PATH_TO_TISSUE_ANALYSIS + 'degreeOutliers.analysis', 'w')
    f.write('%s' % str(degreeOutlierDict))
    f.close()
    f = open(PATH_TO_TISSUE_ANALYSIS + 'vertexOutliers.analysis', 'w')
    f.write('%s' % str(vbOutlierDict))
    f.close()

    return

# Assumes will analyze genemania subgraph with respect to primary_gene_id. 
# Also assumes analyzing PPI network with a reliable name format
def analyzeGenemaniaTissueSubgraph(tissue):
    # Get Network Graph
    inFilename = GENEMANIA_TISSUE_BASE_FILENAME % (tissue, DESIRED_INTERACTION_TYPE)
    graph = readGenemaniaGraph(PATH_TO_TISSUE_SUBGRAPHS + inFilename)
    graph, reportDict = analyzeGraphToolObject(graph)
    return graph, reportDict

def analyzeGraphToolObject(graph):
    # Get V
    numVertices = getNumVertices(graph)

    # Get E
    numEdges = getNumEdges(graph)

    # Get Density
    density = getDensity(graph)

    # Get Connectedness
    connectedness = getConnectedness(graph)

    # Compute Vertex/Edge Betweenness
    computeBetweenness(graph)

    # Get centralDominance
    centralDominance = getCentralPointDominance(graph)

    # Get globalClusteringCoefficient
    globalClusteringCoefficient = getGlobalClusteringCoefficient(graph)

    # Compute localClusteringCoefficients
    computeLocalClusteringCoefficients(graph)

    # Compute degree histogram
    degreeBuckets, degreeCounts = getVertexPropertyHistogram(graph, 'out')

    # Compute Avg Degree
    averageDegree = getDegreeAverage(graph)

    # Get Shortest Distance Histogram
    distanceBuckets, distanceCounts = getShortestDistanceHistogram(graph)

    # Get Diameter
    diameter = getDiameter(graph)

    # Get Bipartiteness
    bipartiteness = getBipartiteness(graph)

    reportDict = {
        'numVertices' : numVertices,
        'numEdges' : numEdges,
        'density' : density,
        'connectedness' : connectedness,
        'centralDominance' : centralDominance,
        'globalClusteringCoefficient' : globalClusteringCoefficient,
        'degreeBuckets' : degreeBuckets,
        'degreeCounts' : degreeCounts,
        'averageDegree' : averageDegree,
        'distanceBuckets' : distanceBuckets,
        'distanceCounts' : distanceCounts,
        'diameter' : diameter,
        'bipartiteness' : bipartiteness
        }

    return graph, reportDict





def scratch():

# Get Betweenness
    vp, ep = graph_tool.centrality.betweenness(graph)

    vpVals = vp.get_array().tolist()

    outfile = open('testout', 'w')
    for i in range(graph.num_vertices()):
        outfile.write('%d\n' % graph.vertex(i).out_degree())
    
    return

    graph_tool.draw.graph_draw(graph, vertex_fill_color=vp,
                          vertex_size=graph_tool.draw.prop_to_size(vp, mi=1, ma=5),
                          vorder=vp, output='cd4_betweenness.png')

#    g = graph_tool.GraphView(g, vfilt)
    sfdp_pos = graph_tool.draw.sfdp_layout(graph)
    frl_pos = graph_tool.draw.fruchterman_reingold_layout(graph, n_iter=1000)
    arf_pos = graph_tool.draw.arf_layout(graph, max_iter=0)

    graph_tool.draw.graphviz_draw(graph, output='graphviz_cd4.png')
#    graph_draw(graph, pos=sfdp_pos, output_size=(1000,1000), output='cd4_sfdp_1000.png')
 #   graph_draw(graph, pos=frl_pos, output_size=(1000,1000), output='cd4_frl_1000.png')
#    graph_draw(graph, pos=arf_pos, output_size=(1000,1000), output='cd4_arf_1000.png')
    graph_draw(graph, pos=sfdp_pos, output='cd4_sfdp.png')
    graph_draw(graph, pos=frl_pos, output='cd4_frl.png')
    graph_draw(graph, pos=arf_pos, output='cd4_arf.png')

# MAIN
def main():
    constructHousekeepingGenesDB()
    fillHousekeepingGenesDB()

#main()
