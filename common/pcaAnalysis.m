% pcaAnalysis.m
% Author: Santhosh Balasubramanian
% Created: April 11, 2013
% Last Modified: April 11, 2013

function [] = pcaAnalysis( inFilePath, inFileName )
    
    % Get Matrix
    M = importdata(strcat(inFilePath, inFileName));
    
    % Perform PCA
    [ D, V ] = pca( M );
    
    % Get Principal Eigenvalues/vectors
    
end