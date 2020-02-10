function [G, D] = bruteslow2(X, k, verbose)
%========================================================================== 
% Syntax
%       [G, D] = bruteslow2(X, k, verbose);
%
% BRUTESLOW2 - Handles the base case of the recursion (Approximate KNN)
%       The implementation of the algorithm is described in: "Fast 
%       Approximate KNN Graph Construction for High Dimensional Data via 
%       Recursive Lanczos Bisection" by Jie Chen, Haw-ren Fang, and Yousef
%       Saad.
%==========================================================================
% Inputs: 
%       X   - data vectors.
%       k   - numbers of nearest neighbors.
%       verbose - (?).
%
% Outputs:   
%       G	- graph.
%       D  	- distances.
%
% Uses L2distances
%==========================================================================
% Reference : "Fast Approximate KNN Graph Construction for High Dimensional
%              Data via Recursive Lanczos Bisection"
% Author   	: Avner Halevy
% Created	: Unknown
% Revised	: Dec 23, 2014 at 09:09 by Karamatou Yacoubou Djima
%==========================================================================

[~,n] = size(X);

if verbose; 
    fprintf('Brute!  : |X|=%d\n',n); 
end

G = zeros(k+1,n);
D = zeros(k+1,n);
labels = X(1,:);

X = X(2:end,:);
dist = L2distance(X,X);

temp = (0:(n-1))*n;

for kk = 1:k+1
    if kk == 1
        II = 1:n;
        M = 1e-10;
    else
        [M, II] = min(dist,[],1);
    end
    G(kk,:) = labels(II);
    D(kk,:) = M;
    
    ind = II + temp;
    dist(ind) = inf;
end

D(1,:) = 1e-10;
D(D==0) = 1e-9;

