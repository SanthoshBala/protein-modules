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
    
    handle = scatter(X, Y, 25, 'filled');
    
    % Marker
    set(handle, 'MarkerEdgeColor', [0.2, 0.2, 1.0]);
    set(handle, 'MarkerFaceColor', [0.2, 0.2, 1.0]);
    
    % Title
    title('Scatter Plot', 'fontname', 'Palatino', 'fontsize', 20);
       
    % X Axis
    xlabel('X Axis', 'fontname', 'Palatino', 'fontsize', 20);
   
    % Y Axis
    ylabel('Y Axis', 'fontname', 'Palatino', 'fontsize', 20);
    
    % Axes
    set(gca, 'fontname', 'Palatino', 'fontsize', 20);
    
    % Line Width
    
    % Line Style
    
    % Axis Label Fonts
    
    % Title
        
    % Log Scale: linear | log
    
end