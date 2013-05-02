% interactionTissueMultiplicity.m
% Author: Santhosh Balasubramanian
% Created: April 16, 2013
% Last Modified: April 16, 2013

PATH_TO_TISSUE_MULTIPLICITY = '/home/santhosh/workspace/thesis/data/analysis/tissue_multiplicity/';
PATH_TO_FIGURES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/';

% Colors
THESIS_RED = [0.8, 0.0, 0.0];
THESIS_ORANGE = [0.92, 0.48, 0.0];
THESIS_YELLOW = [0.9, 0.75, 0.0];
THESIS_LIGHT_GREEN = [0.2, 0.8, 0.2];
THESIS_GREEN = [0.0, 0.49, 0.0];
THESIS_LIGHT_BLUE = [0.02, 0.64, 0.8];
THESIS_BLUE = [0.2, 0.2, 1.0];
THESIS_PURPLE = [0.6, 0.2, 0.6];

% Get Data
M = importdata(strcat(PATH_TO_TISSUE_MULTIPLICITY, ...
        'interaction.tissue.multiplicity'));

handle = figure('position', [100 0 1100 1100], 'paperpositionmode', 'auto', ...
            'color', 'none', 'InvertHardCopy', 'off');
    
% Make Plot
Y = M;
bar(M, 'stack')
    
% Make Colors
patch = findobj(gca, 'Type', 'patch');
colors = [THESIS_PURPLE; THESIS_BLUE; THESIS_LIGHT_BLUE; THESIS_GREEN;...
        THESIS_LIGHT_GREEN; THESIS_YELLOW; THESIS_ORANGE; THESIS_RED];

% Legend    
for i = 1:length(patch)
    set(patch(i), 'FaceColor', colors(i,:));
end
    
% X Axis
xlabel('PPI Tissue Multiplicity', 'fontname', 'Palatino', 'fontsize', 24);
xLabels = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79'};
set(gca, 'xticklabel', xLabels);

% Title
title('Human PPI Tissue Multiplicity', 'fontname', 'Palatino', 'fontsize', 24);

% Y Axis
ylabel('Number of Interactions', 'fontname', 'Palatino', 'fontsize', 24);
set(gca, 'YTickLabel', num2str(transpose(get(gca, 'YTick'))));

% Axes
set(gca, 'fontname', 'Palatino', 'fontsize', 24);
box off;


print(handle, '-depsc2', '-painters', strcat(PATH_TO_FIGURES, 'interaction-tissue-multiplicity.eps'));