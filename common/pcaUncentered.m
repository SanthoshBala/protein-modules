% PCA Uncentered



% Generate ellipse with random data
A = [1 1;
     1.5 2;
     2 1.5;
     2 2;
     2 3;
     3 2;
     3 3;
     3 4;
     4 3;
     4 4;
     4 4.5;
     4.5 4;
     5 5;
     1 3;
     2 4;
     3 5;
     3 1;
     4 2;
     5 3;
    ];
A = A + 0.1*randn(size(A));
X = A(:, 1);
Y = A(:, 2);

handle = figure('position', [100 0 1100 1100], 'paperpositionmode', 'auto', ...
            'color', 'none', 'InvertHardCopy', 'off');
scatter(X, Y, 25, 'filled', 'MarkerEdgeColor', [0.2, 0.2, 1.0], 'MarkerFaceColor', [0.2, 0.2, 1.0]);

% Title
title('Scatter Plot', 'fontname', 'Palatino', 'fontsize', 20);
axis([0 6 0 6]);
% X Axis
xlabel('X Axis', 'fontname', 'Palatino', 'fontsize', 20);

% Y Axis
ylabel('Y Axis', 'fontname', 'Palatino', 'fontsize', 20);

% Axes
set(gca, 'fontname', 'Palatino', 'fontsize', 20);
set(gca, 'xtick', []);
set(gca, 'xticklabel', []);
set(gca, 'ytick', []);
set(gca, 'yticklabel', []);
% Line Width

% Line Style

% Axis Label Fonts

print(handle, '-depsc2', '-painters', '/home/santhosh/plotname.eps');
% Title

% Log Scale: linear | log
