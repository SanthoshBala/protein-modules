#! /usr/bin/python

# networkPCA.py
# Author: Santhosh Balasubramanian
# Created: April 11, 2013
# Last Modified: April 11, 2013


# Python Imports
from copy import *
from random import *

# Library Imports
from pymongo import *

# Global Imports
from settings import *

# Graph Imports
from graphs.graphIO import *
from graphs.graphUtil import *

# Utility Imports


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - #


PATH_TO_PCA_ANALYSIS = PATH_TO_ANALYSIS + 'pca/'
GENE_ASSAY_FILENAME = 'tissue.gene.assay.matrix'
SHUFFLED_GENE_ASSAY_FILENAME = 'shuffled.tissue.gene.assay.matrix'
PPI_ASSAY_FILENAME = 'tissue.ppi.assay.matrix'
SHUFFLED_PPI_ASSAY_FILENAME = 'shuffled.tissue.ppi.assay.matrix'
MODULE_ASSAY_FILENAME = 'tissue.module.assay.matrix'
SHUFFLED_MODULE_ASSAY_FILENAME = 'shuffled.tissue.module.assay.matrix'

# Tissue List with Random Order
RANDOM_TISSUE_LIST = [
    'peripheral_blood_CD8_T_cell',
    'placenta',
    'tonsil',
    'pons',
    'molt4_lymphoblastic_leukemia',
    'olfactory_bulb',
    'colorectal_adenocarcinoma',
    'skin',
    'bone_marrow_CD34',
    'spinal_cord',
    'peripheral_blood_CD19_B_cell',
    'superior_cervical_ganglion',
    'thymus',
    'thyroid',
    'peripheral_blood_CD14_monocyte',
    'temporal_lobe',
    'pancreas',
    'testis_seminiferous_tubule',
    'occipital_lobe',
    'k562_chronic_myelogenous_leukemia',
    'blood',
    'trigeminal_ganglion',
    'dorsal_root_ganglion',
    'bone_marrow_CD33_myeloid',
    'b_lymphoblast',
    'adrenal_cortex',
    'bone_marrow',
    'globus_pallidus',
    'fetal_liver',
    'peripheral_blood_CD4_T_cell',
    'peripheral_blood_CD56_NK_cell',
    'liver',
    'whole_brain',
    'uterus_corpus',
    'prefrontal_cortex',
    'adrenal_gland',
    'cingulate_cortex',
    'brain_amygdala',
    'testis_germ_cell',
    'testis',
    'smooth_muscle',
    'fetal_thyroid',
    'fetal_brain',
    'subthalamic_nucleus',
    'prostate',
    'cerebellar_peduncle',
    'cardiac_myocyte',
    'fetal_lung',
    'testis_interstitial',
    'daudi_burkitts_lymphoma',
    'lung',
    'raji_burkitts_lymphoma',
    'tongue',
    'brain_caudate_nucleus',
    'brain_thalamus',
    'peripheral_blood_BDCA4_dendritic_cell',
    'hl60_promyelocytic_leukemia',
    'bone_marrow_CD71_early_erythroid',
    'skeletal_muscle_psoas',
    'salivary_gland',
    'bronchial_epithelial_cell',
    'medulla_oblongata',
    'parietal_lobe',
    'lymph_node',
    'trachea',
    'heart',
    'islet_cell',
    'hypothalamus',
    'ovary',
    'adipocyte',
    'bone_marrow_CD105_endothelium',
    'uterus',
    'ciliary_ganglion',
    'atrioventricular_node',
    'testis_leydig_cell',
    'appendix',
    'cerebellum',
    'pituitary',
    'kidney'    
    ]


# - - - - - - - - - - ASSAY MATRIX CREATION - - - - - - - - - - #


# createGeneAssayMatrix: Creates matrix describing assays of tissues
# in each gene "condition." For use in PCA analysis.
def createGeneAssayMatrix():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Output File
    outFile = open(PATH_TO_PCA_ANALYSIS + GENE_ASSAY_FILENAME, 'w')

    # Iterate through <ngtmDB>
    for geneRecord in ngtmDB.find():

        geneTissues = geneRecord.get('tissue_list')

        # Iterate through Tissues
        numTissues = len(RANDOM_TISSUE_LIST)
        for i in range(numTissues):
            tissue = RANDOM_TISSUE_LIST[i]

            # Assign Assay Value
            if tissue in geneTissues:
                value = 1
            else:
                value = -1

            # Make Values Tab-Separated
            if i == (numTissues - 1):
                outFile.write('%d' % value)
            else:
                outFile.write('%d\t' % value)

        outFile.write('\n')

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return


# createPPIAssayMatrix: Creates matrix describing assays of tissues
# in each PPI "condition." For use in PCA analysis.
def createPPIAssayMatrix():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Output File
    outFile = open(PATH_TO_PCA_ANALYSIS + GENE_ASSAY_FILENAME, 'w')

    outFile = open(PATH_TO_PCA_ANALYSIS + PPI_ASSAY_FILENAME, 'w')

    # Get Global Graph Interaction Set
    graph = getTissueSubgraph('global')
    interactionSet = getInteractionSet(graph)

    # Iterate through Edges
    for interaction in interactionSet:
        geneA = interaction[0]
        geneB = interaction[1]

        # Get Tissue Lists
        recordA = ngtmDB.find_one( { 'gene_id' : geneA } )
        recordB = ngtmDB.find_one( { 'gene_id' : geneB } )
        if not recordA or not recordB:
            continue
        tissuesA = recordA.get('tissue_list')
        tissuesB = recordB.get('tissue_list')

        # Iterate through Tissues
        numTissues = len(RANDOM_TISSUE_LIST)
        for i in range(numTissues):
            tissue = RANDOM_TISSUE_LIST[i]
            
            # Get Assay Value
            if (tissue in tissuesA) and (tissue in tissuesB):
                value = 1
            else:
                value = -1

            # Make Values Tab-Separated
            if i == (numTissues - 1):
                outFile.write('%d' % value)
            else:
                outFile.write('%d\t' % value)
            
        outFile.write('\n')

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return
        
        
# createModuleAssayMatrix: Creates matrix describing assays of tissues
# in each module "condition." For use in PCA analysis.
def createModuleAssayMatrix():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    modDB = db.modules

    # Open Output File
    outFile = open(PATH_TO_PCA_ANALYSIS + GENE_ASSAY_FILENAME, 'w')

    numModules = modDB.count()

    # Iterate through Modules
    for modRecord in modDB.find():
        modTissues = modRecord.get('tissue_list')

        # Iterate through Tissues
        numTissues = len(RANDOM_TISSUE_LIST)
        for i in range(numTissues):
            tissue = RANDOM_TISSUE_LIST[i]
            
            # Get Assay Value
            if tissue in modTissues:
                value = 1
            else:
                value = -1

            # Make Values Tab-Separated
            if i == (numTissues - 1):
                outFile.write('%d' % value)
            else:
                outFile.write('%d\t' % value)

        outFile.write('\n')

    # Close File and DB Connection
    outFile.close()
    cli.close()

    return


# - - - - - - - - - - SHUFFLED MATRIX CREATION - - - - - - - - - - #



# createShuffledGeneAssayMatrix: Creates gene assay matrix using shuffled
# tissue labels by aligning FUNCTIONAL and RANDOM tissue lists.
def createShuffledGeneAssayMatrix():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    ngtmDB = db.normGeneTissueMap

    # Open Output File
    outFile = open(PATH_TO_PCA_ANALYSIS + SHUFFLED_GENE_ASSAY_FILENAME, 'w')

    # Iterate through <ngtmDB>
    for geneRecord in ngtmDB.find():
        
        geneTissues = geneRecord.get('tissue_list')

        # Iterate through Tissues
        numTissues = len(RANDOM_TISSUE_LIST)
        for i in range(numTissues):
            tissue = Ran
            
        
