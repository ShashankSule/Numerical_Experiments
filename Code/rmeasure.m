% Load images
X = double(rgb2gray(X));

% Parameters setting and space allocation
nscales = 4;
sigmaval = 10;
thresholdingFactor = 1;
SH_I = zeros(size(X,1),size(X,2),nscales+1);

% Create shearlets
shearletSystem = SLgetShearletSystem2D(0,size(X,1),size(X,2),nscales);

% Decomposition
coeffs = SLsheardec2D(X,shearletSystem);

% Thresholding
coeffs = coeffs.*(abs(coeffs) > thresholdingFactor*reshape(repmat(shearletSystem.RMS,[size(X,1)*size(X,2) 1]),[size(X,1),size(X,2),length(shearletSystem.RMS)])*sigmaval);

% Compute sum of coefficients at each scale
nshearlets = size(X,3);

for i = 1:nscales
    for j = 1:nshearlets
        if shearletSystem.shearletIdxs(j,2,:)==i
            SH_I(:,:,i) = coeffs(:,:,j) + SH_I(:,:,i);
        end
    end
    SH_I(:,:,i) = 2^(5*i/4)*SH_I(:,:,i);
end
SH_I(:,:,nscales+1) = coeffs(:,:,end);

% R-measure
Rmax =  max(SH_I,[],3);

