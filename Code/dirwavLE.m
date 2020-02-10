function [MapsX, Xrec, coeffs] = dirwavLE(X, method, nscales, LE_Opts)
%========================================================================== 
% Syntax
%       [MapsX, Xrec] = dirwavLE(X, method, nscales, LE_Opts);
% 
% dirwavLE -- Directional wavelets followed by LE (Laplacian Eigenmaps)
%       This routine has two main steps:
%       1)  A partial directional wavelet reconstruction, which applies the 
%       transform to an image, select a percentage of the largest 
%       coefficients (in modulus), and set the others to zero. Then, it 
%       takes the inverse transform to obtain a partial reconstruction of 
%       the original image. Several methods are currently available:
%       curvelets, shearlets (soon wavelets and contourlets. 
%       2)  The computation of the Laplacian Eigenmaps given the
%       reconstructed image.
%==========================================================================
% Inputs: 
%       X       - data matrix.
%       method:
%               'Shearlets', 'Curvelets', 'Contourlets' and 'Wavelets'  
%       Wav_Opts:
%               nscales - number of scales             
%       LE_Opts:
%               LE_type - 0 is LE on reconstruction, 1 if LE on
%                             coefficients
%               knn     - number of nearest neighbors
%               ndim   	- number of vectors
%               sigma   - scale in heat kernel
%
% Outputs:   
%       MapsX   - Laplacian Eigenimages
%       Xrec    = reconstructed image
%
% Notes:
%   1)  The LE algorithm requires data in more than 1 band. Make sure that 
%       the image is in RBG (or more bands) form.  
%==========================================================================
% Reference : N/A
% Author   	: Karamatou Yacoubou Djima
% Created	: February 27, 2017 at 11:03
% Revised	: April 26, 2017 at 10:20
%==========================================================================

% Declare sizes and initialize matrices and parameter
[m, n, l] = size(X);
Xrec = zeros(m,n,l);

LE_type = LE_Opts(1);
k = LE_Opts(2);
ndim = LE_Opts(3);
sigmaval = LE_Opts(4);
MapsX = zeros(m,n,ndim); 

% Compute directional wavelet coefficients and reconstructed images
disp(['Computing ',lower(method),' coefficients and reconstructed images'])
    
for j = 1:l
    origimg = double(X(:,:,j));
    if strcmp(method,'Wavelets')
        disp('Not setup yet')
        return;
    elseif strcmp(method,'Shearlets')
        [restored_img, coeffs] = shearrec(origimg, nscales);
    elseif strcmp(method,'Curvelets')
        pctg = wav_Opts(1);
        [restored_img] = curverec(origimg, pctg);
    elseif strcmp(method,'Contourlets')
        disp('Not setup yet')
        return;
    else
        disp('Unrecognized directional method')
    end 
    % Xrec(:,:,j) = restored_img/250;
    Xrec(:,:,j) = restored_img;
end

disp('Computing Laplacian eigenmaps...')

% Laplacian Eigenmaps
if LE_type == 0
    [MappedX, ~] = rungenLE(Xrec, k, ndim, sigmaval);
else
    [MappedX, ~] = rungenLE(coeffs, k, ndim, sigmaval);
end

[~, nsz] = size(MappedX);
% Reshape eigenvectors as images
for j = 1:nsz
    MapsX(:,:,j) = reshape(MappedX(:,j),[m,n]);
end

end
