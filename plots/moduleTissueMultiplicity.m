% moduleTissueMultiplicity.m
% Author: Santhosh Balasubramanian
% Created: April 16, 2013
% Last Modified: April 16, 2013

PATH_TO_TISSUE_MULTIPLICITY = '/home/santhosh/workspace/thesis/data/analysis/tissue_multiplicity/';
PATH_TO_FIGURES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/';

% Get Data
A = importdata(strcat(PATH_TO_TISSUE_MULTIPLICITY, ...
        'module.tissue.multiplicity'));
    
handle = figure('position', [100 0 1100 1100], 'paperpositionmode', 'auto', ...
            'color', 'none', 'InvertHardCopy', 'off');
        
hist(A.data, 79, 'BarWidth', 1.0);

% Change Color
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.2, 0.2, 1.0], 'EdgeColor', [1.0, 1.0, 1.0]);

% Title
title('Human Proteome Module Multiplicity', 'fontname', 'Palatino', 'fontsize', 20);

% X Axis
xlabel('Tissue Multiplicity', 'fontname', 'Palatino', 'fontsize', 20);
% Y Axis
ylabel('Number of Modules', 'fontname', 'Palatino', 'fontsize', 20);

% Axes
set(gca, 'fontname', 'Palatino', 'fontsize', 20);
box off;


print(handle, '-depsc2', '-painters', strcat(PATH_TO_FIGURES, 'protein-tissue-multiplicity.eps'));