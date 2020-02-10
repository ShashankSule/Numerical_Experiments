function [Xrec, coeffs, shearletSystem] = shearrec(X, scales)
%========================================================================== 
% Syntax
%       [Xrec, coeffs, shearletSystem] = shearrec(X, scales, thresholdingFactor)
% 
% shearrec -- SHEARlets REConstruction
%       This routine ompute a 2D shearlet system, select a percentage of 
%       the largest coefficients (in modulus), and set the others to zero. 
%       Then, it takes the inverse transform to obtain a partial 
%       reconstruction of the original image.
%==========================================================================
% Inputs: 
%       X       - n x n x 3 or n x n data matrix, n is dyadic
%
% Outputs:   
%       Xrec    - restored image 
%       coeffs  - shearlets coefficients
%==========================================================================
% Reference : ShearLab3D v1.1 Copyright (c) 
% Author   	: Rafael Reisenhofer 
% Created	: 10/11/2014
% Revised	: April 21, 2017 at 7:43 by Karamatou Yacoubou Djima (adapta-
%             tion of 2D [example)
%==========================================================================

% Intitialize data
[~, ~, p] = size(X);
if p == 1
    X = double(X);
elseif p == 3
    X = double(rgb2gray(X));
else
    disp('Image has wrong format. Must be mx3.')
    return;
end

% Parameters
sigma = 30;

%shearLevels = ceil((1:scales)/2);
% shearLevels = [2, 2, 2, 2, 2];
thresholdingFactor = [0 2.5 2.5 2.5 3.8];
% thresholdingFactor = [2.5, 2.5, 2.5, 3.8];
%directionalFilter = modulate2(dfilters('Haar','d'),'c');
%directionalFilter = modulate2(dfilters('cd','d'),'c');

% Create shearlets
shearletSystem = SLgetShearletSystem2D(0,size(X,1),size(X,2),scales);
% shearletSystem = SLgetShearletSystem2D(0,size(X,1),size(X,2),scales,shearLevels,0,directionalFilter);

% Decomposition
coeffs = SLsheardec2D(X,shearletSystem);

% My Thresholding
% coeffs = coeffs.*(abs(coeffs) > thresholdingFactor*reshape(repmat(shearletSystem.RMS,[size(X,1)*size(X,2) 1]),[size(X,1),size(X,2),length(shearletSystem.RMS)])*sigma);

% Scalewise ard thresholding on shearlet coefficients
for j = 1:shearletSystem.nShearlets
    idx = shearletSystem.shearletIdxs(j,:);
    coeffs(:,:,j) = coeffs(:,:,j).*(abs(coeffs(:,:,j)) >= thresholdingFactor(idx(2) + 1)*sigma*shearletSystem.RMS(j));
end

% Reconstruction
Xrec = SLshearrec2D(coeffs,shearletSystem);

