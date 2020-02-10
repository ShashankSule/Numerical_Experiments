function [MapsX, Xrec] = dirwavPCA(X, method, wav_Opts, PCA_Opts)
%========================================================================== 
% Syntax
%       [MapsX] = dirwavPCA(X, method, wav_Opts, PCA_Opts);
% 
% dirwavPCA -- Directional wavelets followed by PCA (principal Component Analysis)
%       This routine has two main steps:
%       1)  A partial directional wavelet reconstruction, which applies the 
%       transform to an image, select a percentage of the largest 
%       coefficients (in modulus), and set the others to zero. Then, it 
%       takes the inverse transform to obtain a partial reconstruction of 
%       the original image. Several methods are currently available:
%       curvelets, shearlets (soon wavelets and contourlets. 
%       2)  The computation of the principal components given the
%       reconstructed image.
%==========================================================================
% Inputs: 
%       X       - data matrix.
%       method:
%               'Shearlets', 'Curvelets', 'Contourlets' and 'Wavelets'  
%       Wav_Opts:
%               'Shearlets' -- scales, thresfac
%               'Curvelets' -- pctg, percentage of curvelets coefficients
%                              considered              
%       PCA_Opts:
%               ndim    - number of dimensions
%
% Outputs:   
%       MapsX   - Laplacian Eigenmaps
%
% Notes:
%   1)  The PCA algorithm requires data in more than 1 band. Make sure that 
%       the image is in RBG (or more bands) form.  
%==========================================================================
% Reference : N/A
% Author   	: Karamatou Yacoubou Djima
% Created	: April 19, 2017 at 10:55
% Revised	: April 19, 2017 at 10:55
%==========================================================================

% Declare sizes and initialize matrices and parameter
[m, n, l] = size(X);
Xrec = zeros(m,n,l);

ndim = PCA_Opts; 

% Compute directional wavelet coefficients and reconstructed images
disp(strcat({'Computing'},{' '},method,{' '},{'coefficients'},{' '},... 
        {'reconstructed images...'}))
for j = 1:l
    origimg = double(X(:,:,j));
    if strcmp(method,'Wavelets')
        disp('Not setup yet')
    elseif strcmp(method,'Shearlets')
        scales = wav_Opts(1);
        thresfac = wav_Opts(2);
        [restored_img, ~] = shearletsrec(origimg, scales, thresfac);
    elseif strcmp(method,'Curvelets')
        pctg = wav_Opts(1);
        restored_img = curveletsrec(origimg, pctg);
    elseif strcmp(method,'Contourlets')
        disp('Not setup yet')
    else
        disp('Unrecognized directional method')
    end 
    Xrec(:,:,j) = restored_img/250;
end

disp('Computing Principal Components Analysis...')

% Principal Components Analysis
vecXrec = reshape(Xrec,[m*n 3]);
[MappedX, ~] = pca(vecXrec, ndim);
MapsX = reshape(MappedX, [m n ndim]);

% Signal end of code    
load gong.mat; sound(y)
disp('Done!');

