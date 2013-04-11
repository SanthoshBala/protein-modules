#! /usr/bin/python

# settings.py
# Author: Santhosh Balasubramanian
# Created: January 16, 2013
# Last Modified: March 24, 2013

# - - - - - - - - - - MACHINE CONFIGURATION - - - - - - - - - #




# - - - - - - - - - - RUN CONFIGURATION - - - - - - - - - #

PRINT_PROGRESS = True
MODULE_SIZE_THRESHOLD = 4

# - - - - - - - - - - FILE PATHS - - - - - - - - - - - #

PATH_TO_THESIS = '/home/santhosh/Dropbox/Thesis/'

PATH_TO_SCRATCH_SPACE = '/home/santhosh/workspace/thesis/'

PATH_TO_CODE = PATH_TO_THESIS + 'code/'

PATH_TO_COMMON_CODE = PATH_TO_CODE + 'common/'
PATH_TO_LIBRARY_CODE = PATH_TO_CODE + 'lib/'

PATH_TO_DATA = PATH_TO_SCRATCH_SPACE + 'data/'

PATH_TO_GRAPHS_CODE = PATH_TO_CODE + 'graphs/'

PATH_TO_SPICI = PATH_TO_LIBRARY_CODE + 'SPICi/'
PATH_TO_SPICI_BINARY = PATH_TO_SPICI + 'src/./spici'



PATH_TO_ANALYSIS = PATH_TO_DATA + 'analysis/'
PATH_TO_PAIRWISE_ANALYSIS = PATH_TO_ANALYSIS + 'pairwise/'
PATH_TO_GLOBAL_LOCAL_DIFF = PATH_TO_ANALYSIS + 'global_local_diff/'
PATH_TO_TISSUE_ANALYSIS = PATH_TO_ANALYSIS + 'tissues/'
PATH_TO_MODULE_ANALYSIS = PATH_TO_ANALYSIS + 'modules/'
PATH_TO_MODULE_TOPOLOGIES = PATH_TO_ANALYSIS + 'module_topology/'
PATH_TO_MODULE_SIMILARITY = PATH_TO_ANALYSIS + 'module_similarity/'
PATH_TO_EXPRESSION_ANALYSIS = PATH_TO_ANALYSIS + 'expression/'


PATH_TO_EXPRESSION = PATH_TO_DATA + 'expression/'
PATH_TO_NETWORKS = PATH_TO_DATA + 'networks/'
PATH_TO_GLOBAL_GRAPH = PATH_TO_NETWORKS + 'genemania/Homo_sapiens.COMBINED/'
PATH_TO_TISSUE_SUBGRAPHS = PATH_TO_GLOBAL_GRAPH + 'tissues/'
PATH_TO_COLLAPSED_TISSUE_SUBGRAPHS = PATH_TO_GLOBAL_GRAPH + 'collapsed_tissues/'
PATH_TO_CANONICAL_TISSUE_SUBGRAPHS = PATH_TO_TISSUE_SUBGRAPHS

PATH_TO_MODULES = PATH_TO_DATA + 'modules/'
PATH_TO_TISSUE_SPICI_MODULES = PATH_TO_MODULES + 'tissue_spici_modules/'
PATH_TO_TISSUE_MODULE_IDS = PATH_TO_MODULES + 'tissue_module_ids/'
PATH_TO_NOMENCLATURE = PATH_TO_DATA + 'nomenclature/'
PATH_TO_ANNOTATION = PATH_TO_DATA + 'annotation/'

PATH_TO_TOPOLOGY_ANALYSIS = PATH_TO_ANALYSIS + 'module_topology/'
PATH_TO_ONTOLOGY_EXPANSION = PATH_TO_TOPOLOGY_ANALYSIS + 'ontology_expansion/'
PATH_TO_FUNCTION_EXPANSION = PATH_TO_ONTOLOGY_EXPANSION + 'functions/'
PATH_TO_PROCESS_EXPANSION = PATH_TO_ONTOLOGY_EXPANSION + 'processes/'

# - - - - - - - - - - FILE NAMES - - - - - - - - - - - - #
EXECUTE_OUTPUT_FILENAME = 'execute.thesis.output'
EXECUTE_ERROR_FILENAME = 'execute.thesis.error'





DESIRED_ORGANISM = 'Homo sapiens'



GENEMANIA_GLOBAL_BASE_FILENAME = 'combined.%s.network'
GENEMANIA_TISSUE_BASE_FILENAME = '%s.combined.%s.network'
GENEMANIA_COLLAPSED_TISSUE_BASE_FILENAME = '%s.combined.%s.collapsed.network'
CANON_GENEMANIA_TISSUE_BASE_FILENAME = 'canon.%s.combined.%s.network'

TISSUE_SPICI_MODULES_BASE_FILENAME = '%s.modules'
TISSUE_MODULE_IDS_BASE_FILENAME = '%s.module.ids'

MODULE_TOPOLOGY_BASE_FILENAME = '%s.module.topology'
FUNCTION_EXPANSION_BASE_FILENAME = '%s.function.expansion'
PROCESS_EXPANSION_BASE_FILENAME = '%s.process.expansion'
MODULE_SIMILARITY_BASE_FILENAME = 'similarity.%s.module.list'

# - - - - - - - - - - DATABASE MANAGEMENT - - - - - - - - - - #


GENEMANIA_INTERACTION_TYPES = [
    'all',
    'coexpression',
    'colocalization',
    'genetic',
    'pathway',
    'ppi',
    'predicted',
    'protein_domains'
    ]

DESIRED_INTERACTION_TYPE = 'ppi'

DESIRED_NOMENCLATURE = 'primary_gene_id'


# Microarray Analysis Settings

ARRAY_EXPRESSION_THRESHOLD = 200

# - - - - - - - - - - TISSUE NOMENCLATURE - - - - - - - - - - #



CANONICAL_TISSUE_LIST = [
    'adipocyte',
    'adrenal_cortex',
    'adrenal_gland',
    'appendix',
    'atrioventricular_node',
    'b_lymphoblast',
    'blood',
    'bone_marrow',
    'bone_marrow_CD105_endothelium',
    'bone_marrow_CD33_myeloid',
    'bone_marrow_CD34',
    'bone_marrow_CD71_early_erythroid',
    'brain_amygdala',
    'brain_caudate_nucleus',
    'brain_thalamus',
    'bronchial_epithelial_cell',
    'cardiac_myocyte',
    'cerebellar_peduncle',
    'cerebellum',
    'ciliary_ganglion',
    'cingulate_cortex',
    'colorectal_adenocarcinoma',
    'daudi_burkitts_lymphoma',
    'dorsal_root_ganglion',
    'fetal_brain',
    'fetal_liver',
    'fetal_lung',
    'fetal_thyroid',
    'globus_pallidus',
    'heart',
    'hl60_promyelocytic_leukemia',
    'hypothalamus',
    'islet_cell',
    'k562_chronic_myelogenous_leukemia',
    'kidney',
    'liver',
    'lung',
    'lymph_node',
    'medulla_oblongata',
    'molt4_lymphoblastic_leukemia',
    'occipital_lobe',
    'olfactory_bulb',
    'ovary',
    'pancreas',
    'parietal_lobe',
    'peripheral_blood_BDCA4_dendritic_cell',
    'peripheral_blood_CD14_monocyte',
    'peripheral_blood_CD19_B_cell',
    'peripheral_blood_CD4_T_cell',
    'peripheral_blood_CD56_NK_cell',
    'peripheral_blood_CD8_T_cell',
    'pituitary',
    'placenta',
    'pons',
    'prefrontal_cortex',
    'prostate',
    'raji_burkitts_lymphoma',
    'salivary_gland',
    'skeletal_muscle_psoas',
    'skin',
    'smooth_muscle',
    'spinal_cord',
    'subthalamic_nucleus',
    'superior_cervical_ganglion',
    'temporal_lobe',
    'testis',
    'testis_germ_cell',
    'testis_interstitial',
    'testis_leydig_cell',
    'testis_seminiferous_tubule',
    'thymus',
    'thyroid',
    'tongue',
    'tonsil',
    'trachea',
    'trigeminal_ganglion',
    'uterus',
    'uterus_corpus',
    'whole_brain'
    ]

FUNCTIONAL_TISSUE_LIST = [
    'cardiac_myocyte',
    'atrioventricular_node',
    'heart',
    'blood',
    'molt4_lymphoblastic_leukemia',
    'hl60_promyelocytic_leukemia',
    'k562_chronic_myelogenous_leukemia',
    'fetal_liver',
    'liver',
    'tongue',
    'salivary_gland',
    'appendix',
    'colorectal_adenocarcinoma',
    'kidney',
    'islet_cell',
    'pancreas',
    'pituitary',
    'adrenal_cortex',
    'adrenal_gland',
    'thyroid',
    'fetal_thyroid',
    'adipocyte',
    'peripheral_blood_CD4_T_cell',
    'peripheral_blood_CD8_T_cell',
    'peripheral_blood_CD14_monocyte',
    'peripheral_blood_CD19_B_cell',
    'peripheral_blood_CD56_NK_cell',
    'peripheral_blood_BDCA4_dendritic_cell',
    'thymus',
    'skin',
    'bone_marrow_CD33_myeloid',
    'bone_marrow_CD34',
    'bone_marrow_CD71_early_erythroid',
    'bone_marrow_CD105_endothelium',
    'bone_marrow',
    'b_lymphoblast',
    'lymph_node',
    'daudi_burkitts_lymphoma',
    'raji_burkitts_lymphoma',
    'tonsil',
    'smooth_muscle',
    'skeletal_muscle_psoas',
    'spinal_cord',
    'pons',
    'brain_caudate_nucleus',
    'globus_pallidus',
    'dorsal_root_ganglion',
    'trigeminal_ganglion',
    'ciliary_ganglion',
    'superior_cervical_ganglion',
    'cingulate_cortex',
    'prefrontal_cortex',
    'medulla_oblongata',
    'brain_amygdala',
    'brain_thalamus',
    'subthalamic_nucleus',
    'hypothalamus',
    'temporal_lobe',
    'occipital_lobe',
    'parietal_lobe',
    'cerebellar_peduncle',
    'cerebellum',
    'whole_brain',
    'fetal_brain',
    'trachea',
    'bronchial_epithelial_cell',
    'lung',
    'fetal_lung',
    'olfactory_bulb',
    'ovary',
    'placenta', 
    'uterus',
    'uterus_corpus',
    'testis_germ_cell',
    'testis_leydig_cell',
    'testis_interstitial',
    'testis_seminiferous_tubule',
    'testis',
    'prostate',        
]

# Some tasks may operate on global interaction network in addition to tissue
# subgraphs. In these cases, 'global' is just another tissue.
AUGMENTED_TISSUE_LIST = ['global', 'intersection'] + FUNCTIONAL_TISSUE_LIST


TISSUE_SYSTEM_MAP = {
    'cardiac_myocyte' : 'circulatory',
    'atrioventricular_node' : 'circulatory',
    'heart' : 'circulatory',
    'blood' : 'circulatory',
    'molt4_lymphoblastic_leukemia' : 'circulatory',
    'hl60_promyelocytic_leukemia' : 'circulatory',
    'k562_chronic_myelogenous_leukemia' : 'circulatory',

    'fetal_liver' : 'digestive',
    'liver' : 'digestive',
    'tongue' : 'digestive',
    'salivary_gland' : 'digestive',
    'appendix' : 'digestive',
    'colorectal_adenocarcinoma' : 'digestive',

    'kidney' : 'urinary',

    'islet_cell' : 'endocrine',
    'pancreas' : 'endocrine',
    'pituitary' : 'endocrine',
    'adrenal_cortex' : 'endocrine',
    'adrenal_gland' : 'endocrine',
    'thyroid' : 'endocrine',
    'fetal_thyroid' : 'endocrine',
    'adipocyte' : 'endocrine',

    'peripheral_blood_CD4_T_cell' : 'immune',
    'peripheral_blood_CD8_T_cell' : 'immune',
    'peripheral_blood_CD14_monocyte' : 'immune',
    'peripheral_blood_CD19_B_cell' : 'immune',
    'peripheral_blood_CD56_NK_cell' : 'immune',
    'peripheral_blood_BDCA4_dendritic_cell' : 'immune',
    'thymus' : 'immune',

    'skin' : 'integumentary',

    'bone_marrow_CD33_myeloid' : 'lymphatic',
    'bone_marrow_CD34' : 'lymphatic',
    'bone_marrow_CD71_early_erythroid' : 'lymphatic',
    'bone_marrow_CD105_endothelium' : 'lymphatic',
    'bone_marrow' : 'lymphatic',
    'b_lymphoblast' : 'lymphatic',
    'lymph_node' : 'lymphatic',
    'daudi_burkitts_lymphoma' : 'lymphatic',
    'raji_burkitts_lymphoma' : 'lymphatic',
    'tonsil' : 'lymphatic',

    'smooth_muscle' : 'musculoskeletal',
    'skeletal_muscle_psoas' : 'musculoskeletal',

    'spinal_cord' : 'nervous',
    'pons' : 'nervous',
    'brain_caudate_nucleus' : 'nervous',
    'globus_pallidus' : 'nervous',
    'dorsal_root_ganglion' : 'nervous',
    'trigeminal_ganglion' : 'nervous',
    'ciliary_ganglion' : 'nervous',
    'superior_cervical_ganglion' : 'nervous',
    'cingulate_cortex' : 'nervous',
    'prefrontal_cortex' : 'nervous',
    'medulla_oblongata' : 'nervous',
    'brain_amygdala' : 'nervous',
    'brain_thalamus' : 'nervous',
    'subthalamic_nucleus' : 'nervous',
    'hypothalamus' : 'nervous',
    'temporal_lobe' : 'nervous',
    'occipital_lobe' : 'nervous',
    'parietal_lobe' : 'nervous',
    'cerebellar_peduncle' : 'nervous',
    'cerebellum' : 'nervous',
    'whole_brain' : 'nervous',
    'fetal_brain' : 'nervous',

    'trachea' : 'respiratory',
    'bronchial_epithelial_cell' : 'respiratory',
    'lung' : 'respiratory',
    'fetal_lung' : 'respiratory',
    'olfactory_bulb' : 'respiratory',

    'ovary' : 'reproductive',
    'placenta' : 'reproductive', 
    'uterus' : 'reproductive',
    'uterus_corpus' : 'reproductive',
    'testis_germ_cell' : 'reproductive',
    'testis_leydig_cell' : 'reproductive',
    'testis_interstitial' : 'reproductive',
    'testis_seminiferous_tubule' : 'reproductive',
    'testis' : 'reproductive',
    'prostate' : 'reproductive',
}

CANONICAL_GERM_LAYER_LIST = [
    'endoderm', 
    'mesoderm',
    'ectoderm'
    ]

TISSUE_GERM_LAYER_MAP = {
    'cardiac_myocyte' : 'mesoderm',
    'atrioventricular_node' : 'mesoderm',
    'heart' : 'mesoderm',
    'kidney' : 'mesoderm',
    'adrenal_cortex' : 'mesoderm',
    'adrenal_gland' : 'mesoderm',
    'ovary' : 'mesoderm',
    'placenta' : 'mesoderm', 
    'uterus' : 'mesoderm',
    'uterus_corpus' : 'mesoderm',
    'testis_germ_cell' : 'mesoderm',
    'testis_leydig_cell' : 'mesoderm',
    'testis_interstitial' : 'mesoderm',
    'testis_seminiferous_tubule' : 'mesoderm',
    'testis' : 'mesoderm',
    'prostate' : 'mesoderm',
    'smooth_muscle' : 'mesoderm',
    'skeletal_muscle_psoas' : 'mesoderm',
    'blood' : 'mesoderm',
    'peripheral_blood_CD4_T_cell' : 'mesoderm',
    'peripheral_blood_CD8_T_cell' : 'mesoderm',
    'peripheral_blood_CD14_monocyte' : 'mesoderm',
    'peripheral_blood_CD19_B_cell' : 'mesoderm',
    'peripheral_blood_CD56_NK_cell' : 'mesoderm',
    'peripheral_blood_BDCA4_dendritic_cell' : 'mesoderm',
    'molt4_lymphoblastic_leukemia' : 'mesoderm',
    'hl60_promyelocytic_leukemia' : 'mesoderm',
    'k562_chronic_myelogenous_leukemia' : 'mesoderm',
    'bone_marrow_CD33_myeloid' : 'mesoderm',
    'bone_marrow_CD34' : 'mesoderm',
    'bone_marrow_CD71_early_erythroid' : 'mesoderm',
    'bone_marrow_CD105_endothelium' : 'mesoderm',
    'bone_marrow' : 'mesoderm',
    'adipocyte' : 'mesoderm',
    'b_lymphoblast' : 'mesoderm',
    'lymph_node' : 'mesoderm',
    'daudi_burkitts_lymphoma' : 'mesoderm',
    'raji_burkitts_lymphoma' : 'mesoderm',
    'appendix' : 'mesoderm',



    'fetal_liver' : 'endoderm',
    'liver' : 'endoderm',
    'islet_cell' : 'endoderm',
    'pancreas' : 'endoderm',
    'bronchial_epithelial_cell' : 'endoderm',
    'lung' : 'endoderm',
    'fetal_lung' : 'endoderm',
    'trachea' : 'endoderm',
    'thyroid' : 'endoderm',
    'fetal_thyroid' : 'endoderm',
    'thymus' : 'endoderm',
    'colorectal_adenocarcinoma' : 'endoderm',
    'tonsil' : 'endoderm',

    'skin' : 'ectoderm',
    'spinal_cord' : 'ectoderm',
    'pons' : 'ectoderm',
    'brain_caudate_nucleus' : 'ectoderm',
    'globus_pallidus' : 'ectoderm',
    'dorsal_root_ganglion' : 'ectoderm',
    'trigeminal_ganglion' : 'ectoderm',
    'ciliary_ganglion' : 'ectoderm',
    'superior_cervical_ganglion' : 'ectoderm',
    'cingulate_cortex' : 'ectoderm',
    'prefrontal_cortex' : 'ectoderm',
    'medulla_oblongata' : 'ectoderm',
    'brain_amygdala' : 'ectoderm',
    'brain_thalamus' : 'ectoderm',
    'subthalamic_nucleus' : 'ectoderm',
    'hypothalamus' : 'ectoderm',
    'temporal_lobe' : 'ectoderm',
    'occipital_lobe' : 'ectoderm',
    'parietal_lobe' : 'ectoderm',
    'cerebellar_peduncle' : 'ectoderm',
    'cerebellum' : 'ectoderm',
    'whole_brain' : 'ectoderm',
    'fetal_brain' : 'ectoderm',
    'olfactory_bulb' : 'ectoderm',
    'tongue' : 'ectoderm',
    'pituitary' : 'ectoderm',
    'salivary_gland' : 'ectoderm',
}

CANONICAL_ORGAN_SYSTEMS = {
    'circulatory' : [
        'blood',
        'heart',
        'cardiac_myocyte',
        'molt4_lymphoblastic_leukemia',
        'hl60_promyelocytic_leukemia',
        'k562_chronic_myelogenous_leukemia',
        'atrioventricular_node',
         ],
    'digestive' : [
        'fetal_liver',
        'tongue',
        'liver',
        'salivary_gland',
        'colorectal_adenocarcinoma',
        'appendix',
        ],
    'endocannabinoid' : [
        ],
    'endocrine' : [
        'islet_cell',
        'pancreas',
        'pituitary',
        'thyroid',
        'fetal_thyroid',
        'adrenal_cortex',
        'adrenal_gland',
        'adipocyte',        
        ],
    'immune' : [
        'peripheral_blood_BDCA4_dendritic_cell',
        'peripheral_blood_CD14_monocyte',
        'peripheral_blood_CD19_B_cell',
        'peripheral_blood_CD4_T_cell',
        'peripheral_blood_CD56_NK_cell',
        'peripheral_blood_CD8_T_cell',
        'thymus',
        ],
    'integumentary' : [
        'skin',
        ],
    'lymphatic' : [
        'bone_marrow',
        'bone_marrow_CD105_endothelium',
        'bone_marrow_CD33_myeloid',
        'bone_marrow_CD34',
        'bone_marrow_CD71_early_erythroid',
        'daudi_burkitts_lymphoma',
        'raji_burkitts_lymphoma',
        'lymph_node',
        'b_lymphoblast',
        'tonsil',
        ],
    'musculoskeletal' : [
        'skeletal_muscle_psoas',
        'smooth_muscle',
        ],
    'nervous' : [
        'whole_brain',
        'prefrontal_cortex',
        'trigeminal_ganglion',
        'hypothalamus',
        'spinal_cord',
        'subthalamic_nucleus',
        'superior_cervical_ganglion',
        'temporal_lobe',
        'cerebellar_peduncle',
        'cerebellum',
        'ciliary_ganglion',
        'cingulate_cortex',
        'dorsal_root_ganglion',
        'fetal_brain',
        'brain_amygdala',
        'brain_caudate_nucleus',
        'brain_thalamus',
        'pons',
        'globus_pallidus',
        'occipital_lobe',
        'parietal_lobe',
        'medulla_oblongata',
        ],
    'reproductive' : [
        'uterus',
        'uterus_corpus',
        'testis',
        'testis_germ_cell',
        'testis_interstitial',
        'testis_leydig_cell',
        'testis_seminiferous_tubule',
        'ovary',
        'prostate',        
        'placenta', 
       ],
    'respiratory' : [
        'trachea',
        'lung',
        'fetal_lung',
        'bronchial_epithelial_cell',
        'olfactory_bulb',
        ],
    'urinary' : [
        'kidney',
        ],
    'vestibular' : [
        ]
}

CANONICAL_ORGAN_SYSTEM_LIST = [
    'nervous',
    'respiratory',
    'circulatory',
    'reproductive',    
    'musculoskeletal',
    'integumentary',
    'lymphatic',
    'immune',
    'endocrine',
    'digestive',
    'urinary',
    ]

TISSUE_BLACKLIST = [
    'molt4_lymphoblastic_leukemia',
    'hl60_promyelocytic_leukemia',
    'k562_chronic_myelogenous_leukemia',
    'fetal_liver',
    'fetal_thyroid',
    'daudi_burkitts_lymphoma',
    'raji_burkitts_lymphoma',
    'fetal_brain',
    'placenta', 
    'fetal_lung',
]

TISSUE_GROUPINGS = {
    'primary' : [
        'blood',
        'heart',
        'tongue',
        'liver',
        'appendix',
        'pancreas',
        'pituitary',
        'thyroid',
        'thymus',
        'skin',
        'lymph_node',
        'tonsil',
        'smooth_muscle',
        'whole_brain',
        'prefrontal_cortex',
        'hypothalamus',
        'spinal_cord',
        'temporal_lobe',
        'cerebellum',
        'brain_amygdala',
        'occipital_lobe',
        'parietal_lobe',
        'medulla_oblongata',
        'kidney',
        'uterus',
        'testis',
        'ovary',
        'prostate',        
        'trachea',
        'lung',
        ]
    }


GERM_LAYER_TISSUE_LIST = [
    'thyroid',
    'fetal_thyroid',
    'thymus',
    'tonsil',
    'islet_cell',
    'pancreas',
    'fetal_liver',
    'liver',
    'trachea',
    'bronchial_epithelial_cell',
    'lung',
    'fetal_lung',
    'colorectal_adenocarcinoma',


    'cardiac_myocyte',
    'atrioventricular_node',
    'heart',
    'blood',
    'peripheral_blood_CD4_T_cell',
    'peripheral_blood_CD8_T_cell',
    'peripheral_blood_CD14_monocyte',
    'peripheral_blood_CD19_B_cell',
    'peripheral_blood_CD56_NK_cell',
    'peripheral_blood_BDCA4_dendritic_cell',
    'molt4_lymphoblastic_leukemia',
    'hl60_promyelocytic_leukemia',
    'k562_chronic_myelogenous_leukemia',
    'bone_marrow_CD33_myeloid',
    'bone_marrow_CD34',
    'bone_marrow_CD71_early_erythroid',
    'bone_marrow_CD105_endothelium',
    'bone_marrow',
    'adipocyte',
    'b_lymphoblast',
    'lymph_node',
    'daudi_burkitts_lymphoma',
    'raji_burkitts_lymphoma',
    'appendix',
    'adrenal_cortex',
    'adrenal_gland',
    'kidney',
    'ovary',
    'placenta', 
    'uterus',
    'uterus_corpus',
    'testis_germ_cell',
    'testis_leydig_cell',
    'testis_interstitial',
    'testis_seminiferous_tubule',
    'testis',
    'prostate',
    'smooth_muscle',
    'skeletal_muscle_psoas',



    'olfactory_bulb',
    'tongue',
    'pituitary',
    'salivary_gland',
    'skin',
    'spinal_cord',
    'pons',
    'brain_caudate_nucleus',
    'globus_pallidus',
    'cingulate_cortex',
    'prefrontal_cortex',
    'medulla_oblongata',
    'brain_amygdala',
    'brain_thalamus',
    'subthalamic_nucleus',
    'hypothalamus',
    'temporal_lobe',
    'occipital_lobe',
    'parietal_lobe',
    'cerebellar_peduncle',
    'cerebellum',
    'whole_brain',
    'fetal_brain',
    'dorsal_root_ganglion',
    'trigeminal_ganglion',
    'ciliary_ganglion',
    'superior_cervical_ganglion',
]

