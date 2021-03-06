% pca.m
% Author: Santhosh Balasubramanian
% Created: April 11, 2013
% Last Modified: April 11, 2013

function [ eigenvalues, eigenvectors ] = pca( M )

    % Get Size of Matrix
    [genes, ~] = size(M);
    
    % Mean-Center Matrix
    for i = 1:genes
        M(i,:) = M(i,:) - mean(M(i,:));
    end
    
    % Get Covariance
    C = cov(M);
    
    % Get Eigenvalues/EigenVectors
    [V, D] = eig(C);
    
    eigenvalues = diag(D);
    eigenvectors = V;
    
end
    