% scatterPlot.m
% Author: Santhosh Balasubramanian
% Created: April 12, 2013
% Last Modified: April 12, 2013

function [] = scatterPlot( M )
    [~, columns] = size(M);
    if columns > 3
        error('Bad Data')
    end
    
    X = M(:,1);
    Y = M(:,2);
    
    handle = scatter(X, Y);
    
    % Marker Color
    set(handle, 'MarkerEdgeColor', 'b');
    set(handle, 'MarkerFaceColor', 'b');
    
    % Marker Size;
    set(handle, 'MarkerSize', 10);
    
    % Line Width
    set(handle, 'LineWidth', 2);
    
    % Line Style
    set(handle, 'LineStyle', '.');
    
    % Axis Label Fonts
    set(handle, 'FontName', 'Helvetica');
    
    % Title
    set(handle, 'Title', 'Scatter Plot');
    
    % Log Scale: linear | log
    set(handle, 'XScale', 'linear');
    set(handle, 'YScale', 'linear');
    set(handle, 'ZScale', 'linear');
    
end