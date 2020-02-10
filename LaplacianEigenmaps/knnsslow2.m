function [G, D]= knnsslow2(X, k, alpha, b, verbose)
%========================================================================== 
% Syntax
%       [G, D]= knnsslow2(X, k, alpha, b, verbose);
%
% KNNSLOW2 - Approximate KNN
%       The implementation of the algorithm is described in: "Fast 
%       Approximate KNN Graph Construction for High Dimensional Data via 
%       Recursive Lanczos Bisection" by Jie Chen, Haw-ren Fang, and Yousef
%       Saad.
%==========================================================================
% Inputs: 
%       X   - data vectors.
%       k   - numbers of nearest neighbors.
%       alpha, b, verbose - (option?).
%
% Outputs:   
%       G	- graph
%       D  	- distances
%
% Uses L2distances
%==========================================================================
% Reference : "Fast Approximate KNN Graph Construction for High Dimensional
%              Data via Recursive Lanczos Bisection"
% Author   	: Avner Halevy
% Created	: October 5, 2009
% Revised	: Dec 23, 2014 at 09:37 by Karamatou Yacoubou Djima
%==========================================================================

[~,n] = size(X);

if n < b
    [G, D] = bruteslow2(X,k,verbose);
else
    [X1,X2] = divides(X, alpha, verbose);
    [G1, D1] = knnsslow2(X1, k, alpha, b, verbose);
    [G2, D2] = knnsslow2(X2, k, alpha, b, verbose);
    [G, D] = conquerslow2(G1, G2, D1, D2, k, verbose);
end


