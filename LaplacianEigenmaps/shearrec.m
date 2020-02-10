function [Xrec, coeffs, shearletSystem] = shearrec(X, scales, thresholdingFactor)
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
sigma = 1;

% shearLevels = [2, 2, 2, 2, 2];
%directionalFilter = modulate2(dfilters('Haar','d'),'c');
%directionalFilter = modulate2(dfilters('cd','d'),'c');

% Create shearlets
shearletSystem = SLgetShearletSystem2D(0,size(X,1),size(X,2),scales);
%shearletSystem = SLgetShearletSystem2D(0,size(X,1),size(X,2),scales,shearLevels,0,directionalFilter);

% Decomposition
coeffs = SLsheardec2D(X,shearletSystem);

% Thresholding
coeffs = coeffs.*(abs(coeffs) > thresholdingFactor*reshape(repmat(shearletSystem.RMS,[size(X,1)*size(X,2) 1]),[size(X,1),size(X,2),length(shearletSystem.RMS)])*sigma);

%coeffssort = sort(abs(coeffs),'desc');
%num = numel(coeffs);
%thr = coeffssort(floor(num*thresholdingFactor));
%coeffs = coeffs.*(abs(coeffs)>=thr);

% Reconstruction
Xrec = SLshearrec2D(coeffs,shearletSystem);

