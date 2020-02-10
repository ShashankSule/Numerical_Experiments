function [D, G] = annfusion(X, sz, k, alpha, brute, gamma, verbose, sym)
%========================================================================== 
% Syntax
%       [D, G] = annfusion(X, sz, k, alpha, brute, gamma, verbose, sym);
%
% ANNFUSION - ANNISOTROPIC DIFFUSION (?)
%==========================================================================
% Inputs: 
%       X       - data vectors
%       k       - numbers of nearest neighbors
%       alpha, brute, verbose, sym : options
%
% Output:   
%       D       - distances
%       G       - graph
%
% See kknslow2
%==========================================================================
% Reference : N/A
% Author   	: Unknown
% Created	: Unknown
% Revised	: Dec 23, 2014 at 09:09 by Karamatou Yacoubou Djima
%==========================================================================

if ~exist('verbose','var')
    verbose = false;
end

if ~exist('sym','var')
    sym=true;
end

% X = X';
[~,n] = size(X); % m:spectral dim, n:no. of pixels
labels = 1:n;

if gamma>0
    [x0, y0, z0] = ind2sub([sz(1) sz(2) sz(3)],1:n);
    if gamma < inf
        Y = [labels; X; gamma * x0; gamma * y0; gamma * z0];
    else
        Y = [labels; x0; y0; z0];
    end
else 
    Y = [labels;X];   
end

disp('Finding edges...')
[G, Dist] = knnsslow2(Y, k, alpha, brute, verbose);

disp('Building sparse matrix...')
ix = repmat(1:n, k+1,1);
D = sparse(ix(:), G(:), Dist(:), n, n);

if sym
D = max(D, D');
end

end
