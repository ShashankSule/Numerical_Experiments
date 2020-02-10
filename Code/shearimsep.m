function [imgCurvesRec, imgPointsRec] = shearimsep(varargin)
%========================================================================== 
% Syntax
%       [imgCurvesRec, imgPointsRec] = shearimsep(X, opt);
% 
% shearimsep -- SHEARlets IMage SEParation
%==========================================================================
% Inputs: 
%       X       - data matrix.
%       opt:    shearLevels type at scale 4      
%               if opt = 0: shearlet system sl2d1 (4 scales, redundancy 25)
%               if opt = 1: shearlet system sl2d1 (4 scales, redundancy 49)
%               Default is opt = 1. 
%
% Outputs:   
%       MapsX   - Laplacian Eigenimages
%       Xrec    = reconstructed image
%
% Notes:
%   1)  The LE algorithm requires data in more than 1 band. Make sure that 
%       the image is in RBG (or more bands) form.  
%==========================================================================
% Reference : ShearLab3D v1.1 Copyright (c) 
% Author   	: Shearlab
% Created	: 2014
% Revised	: April 21, 2017 at 8:39 by Karamatou Yacoubou Djima
%==========================================================================


if nargin < 2
    X = varargin{1};
    testSL2D = 1;
else
    X = varargin{1};
    testSL2D = varargin{2};
end  

% Initialize parameters
iterations = 100; % number of iterations during iterative thresholding
gamma = 3; % used for soft thresholding during total variation correction
stopFactor = 3; % determines lowest thresholdg during iterative thresholding
scalesWavelet = 4; % number of scales of the stationary wavelet transform
freqWeights = [.1, .1, 2, 2]; % reiweights the frequency domain using an atrous decomposition
sizeX = size(X,1);
sizeY = size(X,2);

% Determine what type of shearlet system must be used
if testSL2D == 0    % test shearlet system sl2d1 (4 scales, redundancy 25)
    testSL2D1 = 1; 
    testSL2D2 = 0;
else                % test shearlet system sl2d2 (4 scales, redundancy 49)
    testSL2D1 = 0;
    testSL2D2 = 1; 
end

% Make shearlets decomposition
% SL2D1
if testSL2D1
    sl2d1 = SLgetShearletSystem2D(0,sizeX,sizeY,4,[0 0 1 1]);
end
% SL2D2
if testSL2D2
    sl2d2 = SLgetShearletSystem2D(0,sizeX,sizeY,4,[1 1 2 2]);
end

% Make separation
if testSL2D1
    [imgCurvesRec, imgPointsRec] = SLExperimentSeparate(iXg,scalesWavelet,iterations,stopFactor,gamma,freqWeights,sl2d1);

    imgPointsRec = (imgPointsRec - min(imgPointsRec(:)))*(255/max(max(imgPointsRec - min(imgPointsRec(:)))));
    imgCurvesRec = (imgCurvesRec - min(imgCurvesRec(:)))*(255/max(max(imgCurvesRec - min(imgCurvesRec(:)))));        
end
if testSL2D2
    [imgCurvesRec, imgPointsRec] = SLExperimentSeparate(X,scalesWavelet,iterations,stopFactor,gamma,freqWeights,sl2d2);

    imgPointsRec = (imgPointsRec - min(imgPointsRec(:)))*(255/max(max(imgPointsRec - min(imgPointsRec(:)))));
    imgCurvesRec = (imgCurvesRec - min(imgCurvesRec(:)))*(255/max(max(imgCurvesRec - min(imgCurvesRec(:)))));
end
       
