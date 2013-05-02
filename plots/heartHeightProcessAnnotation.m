% heartDegreeAnnotation.m

PATH_TO_ANNOTATION = '/home/santhosh/workspace/thesis/data/analysis/topology/module_topology/';
PATH_TO_FIGURES = '/home/santhosh/Dropbox/Thesis/report/final_report/images/';

% Get Data
data = importdata(strcat(PATH_TO_ANNOTATION, 'heart.module.topology.process.annotation'));

X = data(:,1);
Y = data(:,2);

handle = figure('position', [100 0 1100 1100], 'paperpositionmode', 'auto', ...
            'color', 'none', 'InvertHardCopy', 'off');
scatter(X, Y, 25, 'filled', 'MarkerEdgeColor', [0.2, 0.2, 1.0], 'MarkerFaceColor', [0.2, 0.2, 1.0]);

% Title
%title('Heart Module Topology Process Annotation', 'fontname', 'Palatino', 'fontsize', 20);

% X Axis
xlabel('Vertex Height', 'fontname', 'Palatino', 'fontsize', 20);

% Y Axis
ylabel('# of Process Annotations', 'fontname', 'Palatino', 'fontsize', 20);

% Axes
set(gca, 'fontname', 'Palatino', 'fontsize', 20);
set(gca, 'xticklabel', {'1', '2', '3', '4', '5', '6', '7'});

print(handle, '-depsc2', '-painters', strcat(PATH_TO_FIGURES, 'heart-topology-annotation.eps'));
