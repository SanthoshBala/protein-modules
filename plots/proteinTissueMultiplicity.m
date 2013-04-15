% proteinTissueMultiplicity.m
% Author: Santhosh Balasubramanian
% Created: April 14, 2013
% Last Modified: April 14, 2013


PATH_TO_TISSUE_SPECIFICITY = '/home/santhosh/workspace/thesis/data/analysis/tissue_specificity/';
PATH_TO_FIGURES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/';
% Get Data
A = importdata(strcat(PATH_TO_TISSUE_SPECIFICITY, ...
        'protein.tissue.multiplicity'));
    
handle = figure('position', [100 0 1100 1100], 'paperpositionmode', 'auto', ...
            'color', 'none', 'InvertHardCopy', 'off');
        
hist(A.data, 79);

% Change Color
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.2, 0.2, 1.0], 'EdgeColor', [1.0, 1.0, 1.0]);

% Title
title('Human Proteome Tissue Multiplicity', 'fontname', 'Palatino', 'fontsize', 20);

% X Axis
xlabel('Number of Tissues', 'fontname', 'Palatino', 'fontsize', 20);

% Y Axis
ylabel('Number of Proteins', 'fontname', 'Palatino', 'fontsize', 20);

% Axes
set(gca, 'fontname', 'Palatino', 'fontsize', 20);
box off;


print(handle, '-depsc2', '-painters', strcat(PATH_TO_FIGURES, 'protein-tissue-multiplicity.eps'));