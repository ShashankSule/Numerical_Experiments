function [mappedX, mapping, lambda] = rungenLE(X, k, no_dims, sigmaval)
%========================================================================== 
% Syntax
%       [MappedX, mapping, lambda] = rungenLE(X, k, no_dims, sigmaval);
%
% rungenLE -- RUN GENeral Laplacian Eigenmaps
%       This routine computes the Laplacian Eigenmaps given a data set X.
%==========================================================================
% Inputs: 
%       X       - data matrix
%       k       - number of nearest neighbors
%       no_dims - number of eigenvectors
%       sigma   - heat kernel parameter
%
% Outputs:   
%       MappedX - Eigenmaps
%       Mapping - (?)
%       lambda  - 
% The eigenvectors are computed using JDQR.
%==========================================================================
% Reference : N/A
% Author   	: Unknown
% Created	: Unknown
% Revised	: Dec 23, 2014 at 10:17 by Karamatou Yacoubou Djima
%==========================================================================

gamma = 0;
normalize = 1;

[G, ~] = nngraph(X, k, gamma, normalize);

[mappedX, mapping, lambda] = lapbasic(G, no_dims, sigmaval, 'JDQR');

% Signal end of code    
load gong.mat; sound(y)
disp('Done!');

end
