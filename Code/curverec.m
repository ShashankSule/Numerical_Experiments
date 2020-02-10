function restored_img = curverec(X, pctg)
%========================================================================== 
% Syntax
%       Xrec = curverec(X, pctg)
%
% curverec --  Partial CURVElets REConstruction.
%       This function applies the curvelet transform to an image, select a 
%       percentage of the largest coefficients (in modulus), and set the 
%       others to zero. Then, it takes the inverse curvelet transform to 
%       obtain a partial reconstruction of the original image.
%==========================================================================
% Inputs: 
%       X       - (m x n) image
%       pctg    - percentage of coefficients used in partial reconstruction
%
% Outputs:   
%       Xrec    - (m x n) reconstructed image
%
% Notes: This code is basically the script fdct_wrapping_demo_recon.m made
%       into a function for greater flexibility (parameter change, use in 
%       other functions).
%==========================================================================
% Reference : fdct_usfft_demo_recon.m in Curvelet toolbox
% Author   	: Emmanuel Candès, Laurent Demanet, Lexing Ying
% Created	: California Institute of Technology, 2005-2007
% Revised	: April 05, 2017 at 10:23 by Karamatou Yacoubou Djima
%==========================================================================

cd ~/'Dropbox (Amherst College)'/Research/Computing/LocalCodes/DirectionalWavelets/Curvelets/fdct_usfft_matlab

% Forward curvelet transform
C = fdct_usfft(double(X),0);

% Get threshold value
cfs =[];
for s=1:length(C)
  for w=1:length(C{s})
    cfs = [cfs; abs(C{s}{w}(:))];
  end
end
cfs = sort(cfs); cfs = cfs(end:-1:1);
nb = round(pctg*length(cfs));
cutoff = cfs(nb);

% Set small coefficients to zero
for s=1:length(C)
  for w=1:length(C{s})
    C{s}{w} = C{s}{w} .* (abs(C{s}{w})>cutoff);
  end
end
restored_img = real(ifdct_usfft(C,0));

cd ~/'Dropbox (Amherst College)'/Research/Computing/MATLAB/ProjectBased/VesselEnhancement/Code
