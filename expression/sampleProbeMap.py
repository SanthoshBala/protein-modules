#! /usr/bin/python

# sampleProbeMap.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: March 24, 2013


# Library Imports
from pymongo import *

# Global Imports
from settings import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


PATH_TO_ARRAY_FILE = PATH_TO_EXPRESSION + 'gse1133/geo/'
ARRAY_FILENAME = 'gse1133.array.soft'

# String Constants
DATABASE = 'DATABASE'
PLATFORM = 'PLATfORM'
SERIES = 'SERIES'
SAMPLE = 'SAMPLE'
SAMPLE_TABLE = 'SAMPLE_TABLE'

MIN_SAMPLE_ID = 'GSM18706'

# Array Fields to Canonical Values
CANONICAL_ARRAY_FIELDS = {
    '!Sample_title' : 'source_name',
    '!Sample_geo_accession' : 'geo_id',
    '!Sample_type' : 'type',
    '!Sample_description' : 'original_description',
    '!Sample_organism_ch1' : 'organism',
    '!Sample_taxid_ch1' : 'ncbi_organism_id',
    '!Sample_platform_id' : 'platform_id',
    '!Sample_series_id' : 'series_id',
    '!Sample_label_ch1' : 'label',
    'tissue' : 'tissue',
    'probes': 'probes'
    }

# Maps from Su et al. tissue names to a canonical name used throughout code.
CANONICAL_TISSUE_MAP = {
    'fetalThyroid' : 'fetal_thyroid',
    'FETALTHYROID' : 'fetal_thyroid',
    'fetalbrain' : 'fetal_brain',
    'fetal brain' : 'fetal_brain',
    'fetalliver' : 'fetal_liver',
    'fetal liver' : 'fetal_liver',
    'fetallung' : 'fetal_lung',
    'fetal lung' : 'fetal_lung',
    'Colorectal Adenocarcinoma' : 'colorectal_adenocarcinoma',
    'Colorectal Adenocarc' : 'colorectal_adenocarcinoma',
    'leukemia, chronic Myelogenous K-562':'k562_chronic_myelogenous_leukemia',
    'leukemiachronicmyelogenous(k562)' : 'k562_chronic_myelogenous_leukemia',
    'leukemialymphoblastic(molt4)': 'molt4_lymphoblastic_leukemia',
    'LEUKEMIALYMPHOBLASTIC(MOLT4)': 'molt4_lymphoblastic_leukemia',
    'leukemia lymphoblastic (MOLT-4)': 'molt4_lymphoblastic_leukemia',
    'leukemia, promyelocytic-HL-60': 'hl60_promyelocytic_leukemia',
    'leukemiapromyelocytic(hl60)' : 'hl60_promyelocytic_leukemia',
    'lymphomaburkittsDaudi' : 'daudi_burkitts_lymphoma',
    'LymphomaRaji': 'raji_burkitts_lymphoma',
    '721_BLymphoblasts' : 'b_lymphoblast',
    '721_Blymphoblasts' : 'b_lymphoblast',
    'Heart' : 'heart',
    'HEART' : 'heart',
    'liver' : 'liver',
    'LIVER' : 'liver',
    'Liver' : 'liver',
    'AdrenalCortex' : 'adrenal_cortex',
    'Appendix' : 'appendix',
    'Brain Amygdala': 'brain_amygdala',
    'CardiacMyocytes': 'cardiac_myocyte',
    'CerebellumPeduncles': 'cerebellar_peduncle',
    'CingulateCortex': 'cingulate_cortex',
    'DRG': 'dorsal_root_ganglion',
    'HBEC': 'bronchial_epithelial_cell',
    'HUMANCULTUREDADIPOCYTE': 'adipocyte',
    'Hypothalamus': 'hypothalamus',
    'Islet': 'islet_cell',
    'IsletCell': 'islet_cell',
    'KIDNEY': 'kidney',
    'Lung': 'lung',
    'Medulla_Oblongata': 'medulla_oblongata',
    'OccipitalLobe': 'occipital_lobe',
    'OlfactoryBulb': 'olfactory_bulb',
    'Ovary': 'ovary',
    'PLACENTA': 'placenta',
    'PROSTATE': 'prostate',
    'Pancreas': 'pancreas',
    'ParietalLobe': 'parietal_lobe',
    'Pituitary': 'pituitary',
    'Pons': 'pons',
    'PrefrontalCortex': 'prefrontal_cortex',
    'Prostate': 'prostate',
    'Skeletal_Muscle_Psoas': 'skeletal_muscle_psoas',
    'SmoothMuscle': 'smooth_muscle',
    'Superior_Cervical_Ganglion': 'superior_cervical_ganglion',
    'TONGUE': 'tongue',
    'TemporalLobe': 'temporal_lobe',
    'Testi-GermCell': 'testis_germ_cell',
    'TestiGermCell': 'testis_germ_cell',
    'TestiIntersitial': 'testis_interstitial',
    'TestiLeydigCell': 'testis_leydig_cell',
    'TestiSeminiferousTubule': 'testis_seminiferous_tubule',
    'Testi_GermCell': 'testis_germ_cell',
    'Testi_Intersitial': 'testis_interstitial',
    'Testi_LeydigCell': 'testis_leydig_cell',
    'Testi_SeminiferousTubule': 'testis_seminiferous_tubule',
    'Thyroid': 'thyroid',
    'Tonsil': 'tonsil',
    'Trigeminal_Ganglion': 'trigeminal_ganglion',
    'UTERUS': 'uterus',
    'Uterus': 'uterus',
    'Uterus_Corpus': 'uterus_corpus',
    'WHOLEBLOOD': 'blood',
    'WHOLEBLOOD(JJV)': 'blood',
    'Whole Brain': 'whole_brain',
    'adrenal gland': 'adrenal_gland',
    'adrenalgland': 'adrenal_gland',
    'atrioventricular_node': 'atrioventricular_node',
    'bone marrow': 'bone_marrow',
    'bone marrow-CD105Endothelial': 'bone_marrow_CD105_endothelium',
    'bone marrow-CD105endothelia': 'bone_marrow_CD105_endothelium',
    'bone marrow-CD105endothelial': 'bone_marrow_CD105_endothelium',
    'bone marrow-CD33Myeloid' : 'bone_marrow_CD33_myeloid',
    'bone marrow-CD33myeloid' : 'bone_marrow_CD33_myeloid',
    'bone marrow-CD34' : 'bone_marrow_CD34',
    'bone marrow-CD71EarlyErythroid': 'bone_marrow_CD71_early_erythroid',
    'bonemarrow' : 'bone_marrow',
    'brain, caudate nucleus' : 'brain_caudate_nucleus',
    'brainThalamus' : 'brain_thalamus',
    'caudatenucleus' : 'brain_caudate_nucleus',
    'cerebellum' : 'cerebellum',
    'ciliary_ganglion' : 'ciliary_ganglion',
    'globus_pallidus' : 'globus_pallidus', 
    'kidney' : 'kidney',
    'lung' : 'lung',
    'lymph node' : 'lymph_node',
    'lymphnode' : 'lymph_node',
    'peripheral blood-BDCA4DendriticCells':'peripheral_blood_BDCA4_dendritic_cell',
    'peripheral blood-BDCA4DentriticCells':'peripheral_blood_BDCA4_dendritic_cell',
    'peripheral blood-CD14Monocytes' : 'peripheral_blood_CD14_monocyte',
    'peripheral blood-CD19BCells' : 'peripheral_blood_CD19_B_cell',
    'peripheral blood-CD4TCells' : 'peripheral_blood_CD4_T_cell',
    'peripheral blood-CD56NKCells' : 'peripheral_blood_CD56_NK_cell',
    'peripheral blood-CD8TCells' : 'peripheral_blood_CD8_T_cell',
    'pituitary' : 'pituitary',
    'salivarygland' : 'salivary_gland',
    'skin' : 'skin',
    'spinal cord' : 'spinal_cord',
    'spinalcord' : 'spinal_cord',
    'subthalamic_nucleus' : 'subthalamic_nucleus',
    'testis' : 'testis',
    'thymus' : 'thymus',
    'thyroid' : 'thyroid',
    'trachea' : 'trachea',
    'wholebrain' : 'whole_brain'
    }


# - - - - - - - - - - SAMPLE:PROBE MAP CREATION - - - - - - - - - - #


# createSampleProbeMap: Creates <sampleProbeMap> from microarray.
def createSampleProbeMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    sampleProbeMap = db.sampleProbeMap

    # Open Input File
    inFile = open(PATH_TO_ARRAY_FILE + ARRAY_FILENAME, 'r')

    # Initialize Variables
    currentEntityType = DATABASE
    currentEntityID = ''
    sampleOrganism = ''
    sampleDescription = ''
    sampleTissueID = 0
    numTissues = 0

    skipLine = False

    record = {}
    # Iterate through Array
    for line in inFile:
        # Skip Comments
        if line[0] == '#':
            continue

        # Potentially Skip Line
        if skipLine:
            skipLine = False
            continue

        if line[0] == '^':
            # Parse New Entity
            lineFields = line.split('=')
            lineFields = [ s.strip() for s in lineFields ]

            # Omit the ^ in Entity Type
            currentEntityType = lineFields[0][1:]
            currentEntityID = lineFields[-1]
            continue
        
        # DATABASE, SERIES, and PLATFORM Info Captured by Annotation Files
        if currentEntityType == DATABASE:
            continue
        elif currentEntityType == SERIES:
            continue
        elif currentEntityType == PLATFORM:
            continue
        elif currentEntityType == SAMPLE:
            # If Sample ID not Human, Ignore It
            if currentEntityID < MIN_SAMPLE_ID:
                continue

            if line[0] == '!':
                # Parse Line
                lineFields = line.split('=')
                lineFields = [ s.strip() for s in lineFields ]

                # Check for Beginning of Sample Table
                if lineFields[0] == '!sample_table_begin':
                    currentEntityType = SAMPLE_TABLE
                    record.update( { 'probes' : {} } )
                    # Skip Header in Sample Table
                    skipLine = True

                if lineFields[0] in CANONICAL_ARRAY_FIELDS.keys():
                    key = CANONICAL_ARRAY_FIELDS.get( lineFields[0] )
                    value = lineFields[1]
                    record.update( { key : value } )
                    if lineFields[0] == '!Sample_description':
                        tissue = CANONICAL_TISSUE_MAP.get(value)
                        record.update( { 'tissue' : tissue } )

        elif currentEntityType == SAMPLE_TABLE:
            if '!sample_table_end' in line:
                # Commit Record to DB
                sampleProbeMap.save(record)
                record = {}
                continue

            lineFields = line.split('\t')
            lineFields = [ s.strip() for s in lineFields ]
            record['probes'].update( { lineFields[0] : lineFields[1] } )

    # Close File and DB Connection
    inFile.close()
    cli.close()
    
    return


# - - - - - - - - - - SAMPLE:PROBE MAP NORMALIZATION - - - - - - - - - - #


# normalizeSampleProbeMap: Removes all probes from <sampleProbeMap> that are
# below the expression threshold, yielding <normSampleProbeMap>.
def normalizeSampleProbeMap():
    # Open DB Connection
    cli = MongoClient()
    db = cli.db
    sampleProbeMap = db.sampleProbeMap
    normSampleProbeMap = db.normSampleProbeMap

    # Iterate through <sampleProbeMap>
    for sampleRecord in sampleProbeMap.find():
        probeDict = sample.get('probes')
        deletionList = []
        
        # Iterate through all Probes
        for probe, expression in probeDict.iteritems():
            if float(expression) < ARRAY_EXPRESSION_THRESHOLD:
                deletionList = deletionList + [ probe ]

        # Iterate through Deletion List and Delete
        for probe in deletionList:
            del sampleRecord['probes'][probe]

        # Commit Result to <normSampleProbeMap>
        normSampleProbeMap.save(sampleRecord)

    # Close DB Connection
    cli.close()

    return
