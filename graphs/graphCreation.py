#! /usr/bin/python

# graphCreation.py
# Author: Santhosh Balasubramanian
# Created: April 10, 2013
# Last Modified: April 10, 2013


# Python Imports
import os
import shutil

# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Utility Imports
from common.strings import *


# - - - - - - - - - - MODULE SUBGRAPHS - - - - - - - - - - #


# induceModuleSubgraph
def induceModuleSubgraph(graph, geneList):
    idVertexMap = getIDMapForGraph(graph)

    subgraph = Graph(directed=False)
    subgraph.vertex_properties['gene_id'] = subgraph.new_vertex_property('string')
    geneIDProp = subgraph.vertex_properties['gene_id']

    numGenes = len(geneList)
    # Iterate through first gene in pair
    for i in range(numGenes):
        geneA = geneList[i]
        vertexA = idVertexMap.get(geneA)

        # Iterate through second gene in pair
        for j in range(i, numGenes):
            geneB = geneList[j]
            vertexB = idVertexMap.get(geneB)
            
            edge = graph.edge(vertexA, vertexB)

            if not edge:
                continue
            else:
                subvertexA = getOrCreateVertex(subgraph, geneIDProp, geneA)
                subvertexB = getOrCreateVertex(subgraph, geneIDProp, geneB)
                subgraph.add_edge(subvertexA, subvertexB)

    return subgraph
