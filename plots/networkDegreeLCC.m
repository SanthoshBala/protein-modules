% networkDegreeDistribution.m
% Author: Santhosh Balasubramanian
% Created: April 18, 2013
% Last Modified: April 18, 2013

PATH_TO_DEGREE_LCC = '/home/santhosh/workspace/thesis/data/analysis/topology/lcc/';

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

aVector = zeros(1, length(RANDOM_TISSUE_LIST));
kVector = zeros(1, length(RANDOM_TISSUE_LIST));
rVector = zeros(1, length(RANDOM_TISSUE_LIST));

for i = 1:length(RANDOM_TISSUE_LIST)
    tissue = RANDOM_TISSUE_LIST(i);
    inFileName = strcat(tissue, '.network.degree.lcc');
    data = importdata(char(strcat(PATH_TO_DEGREE_LCC, inFileName)));

    X = data(:,1);
    Y = data(:,2);

    % Remove 0's, which will become -Inf after log
    cleanX = X;
    cleanY = Y;
    removalIndices = [];
    for j = 1:length(Y)
        if Y(j) == 0
            removalIndices = [ removalIndices j ];
        end
    end
    cleanX(removalIndices,:) = [];
    cleanY(removalIndices,:) = [];

    % Get Best Fit
    logX = log(cleanX);
    logY = log(cleanY);
    p = polyfit(logX, logY, 1);
    k = p(1);
    a = exp(p(2));
    
    aVector(i) = a;
    kVector(i) = k;

    % Get R^2
    rStruct = regstats(logX, logY, 'linear', 'rsquare');
    r = rStruct.rsquare;
    
    rVector(i) = r;
    
end

meanA = mean(aVector);
meanK = mean(kVector);
meanR = mean(rVector);

% Write to File
fileID = fopen(strcat(PATH_TO_DEGREE_LCC, 'average.network.degree.lcc'), 'w');
fprintf(fileID, '%s\t%f\n', 'Mean A', meanA);
fprintf(fileID, '%s\t%f\n', 'Mean K', meanK);
fprintf(fileID, '%s\t%f\n', 'Mean R^2', meanR);
fclose(fileID);