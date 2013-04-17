% pcaAnalysis.m
% Author: Santhosh Balasubramanian
% Created: April 11, 2013
% Last Modified: April 11, 2013

function [D, V, RANDOM_TISSUE_LIST] = pcaAnalysis( inFileName )

    PATH_TO_PCA_DATA = '/home/santhosh/workspace/thesis/data/analysis/pca/';

    RANDOM_TISSUE_LIST = {
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
};
    
    % Get Matrix
    M = importdata(strcat(PATH_TO_PCA_DATA, inFileName));
    
    % Perform PCA
    [ D, V ] = pca( M );
    
    % Get Principal Eigenvalues/Vectors
    primaryVector = V(:,79);
    
    for i = 1:length(primaryVector)
        if primaryVector(i) > 0.1
            disp(RANDOM_TISSUE_LIST(i));
        end
    end
end