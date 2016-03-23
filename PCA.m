function [PC, pcVar, pv_explained] = PCA(x)
%PCA performs principal components analysis on raw data
%
%   [PC, pcVar, pv_explained] = PCA(x)
%
%   takes one argument:
%       x - data matrix with observations across columns and conditions
%           down rows
%
%   Erik Schluter - 8/3/2015

    n = size(x,2);
    
    % subtract off mean for each condition (rows)
    x    = x - repmat(mean(x,2),1,n);
    covD = (1/(n-1))*(x*x');
    
    % find eigenvalues and right eigenvectors (PCs)
    [PC, Lambda] = eig(covD);
    
    % extract variances of PCs and sort
    [pcVar, ind] = sort(diag(Lambda), 'descend');
    % sort the PC coeffs
    PC           = PC(:,ind);
    % compute the percentage of variance explained by each PC
    pv_explained = (pcVar/sum(pcVar))*100;
    
return