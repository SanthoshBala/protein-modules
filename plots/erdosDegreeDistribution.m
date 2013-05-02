% heartDegreeDistribution.m
% Author: Santhosh Balasubramanian
% Created: April 18, 2013
% Last Modified: April 18, 2013

PATH_TO_DEGREE_DISTRIBUTION = '/home/santhosh/workspace/thesis/data/analysis/topology/degree/';
PATH_TO_FIGURES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/';

% Get Data
data = importdata(strcat(PATH_TO_DEGREE_DISTRIBUTION, 'heart.network.degree.distribution'));

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


% Get Best Fit -- Poisson
logX = log(cleanX);
logY = log(cleanY);
p = polyfit(logX, logY, 1);
k = p(1)
a = exp(p(2))
plot(X, a*X.^k, 'LineWidth', 2, 'Color', [0.8, 0.0, 0.0]);

hold off;

h = legend({'Data', '$P(k) = 0.78k^{-1.65}$'});
set(h, 'interpreter', 'latex', 'fontname', 'Palatino');
