% heartDegreeLCC.m
% Author: Santhosh Balasubramanian
% Created: April 22, 2013
% Last Modified: April 22, 2013

PATH_TO_DEGREE_DISTRIBUTION = '/home/santhosh/workspace/thesis/data/analysis/topology/lcc/';
PATH_TO_FIGURES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/';

% Get Data
data = importdata(strcat(PATH_TO_DEGREE_DISTRIBUTION, 'heart.network.degree.lcc'));

X = data(:,1);
Y = data(:,2);
cleanX = X;
cleanY = Y;
removalIndices = [];
for i = 1:length(Y)
    if Y(i) == 0
        removalIndices = [ removalIndices i ];
    end
end
cleanX(removalIndices,:) = [];
cleanY(removalIndices,:) = [];


handle = figure('position', [100 0 1100 1100], 'paperpositionmode', 'auto', ...
            'color', 'none', 'InvertHardCopy', 'off');
plot(X, Y, 'o','MarkerSize', 5, 'MarkerEdgeColor', [0.2, 0.2, 1.0], 'MarkerFaceColor', [0.2, 0.2, 1.0]);

% Title
%title('Heart PPI Network LCC vs. Degree (Log-Log Scale)', 'fontname', 'Palatino', 'fontsize', 24);

% X Axis
xlabel('Degree, k', 'fontname', 'Palatino', 'fontsize', 24);
%set(gca, 'yticklabel', {'0.001', '0.01', '0.1', '1.0'});
%set(gca, 'xticklabel', {'1', '10', '100', '1000'});

% Y Axis
ylabel('Local Clustering Coefficient, C(k)', 'fontname', 'Palatino', 'fontsize', 24);

% Axes
set(gca, 'fontname', 'Palatino', 'fontsize', 24);
box off;

hold on;

% Get Best Fit
logX = log(cleanX);
logY = log(cleanY);
p = polyfit(logX, logY, 1);
k = p(1)
a = exp(p(2))
sortData = sortrows([X a*X.^k], [1 2]);
uniqueData = unique(sortData, 'rows')
plot(uniqueData(:,1), uniqueData(:,2), 'LineWidth', 2, 'Color', [0.8, 0.0, 0.0]);


hold off;

h = legend({'Data', '$C(k) = 0.544k^{-0.579}$'});
set(h, 'interpreter', 'latex', 'fontname', 'Palatino');

print(handle, '-depsc2', '-painters', strcat(PATH_TO_FIGURES, 'heart-degree-lcc.eps'));
