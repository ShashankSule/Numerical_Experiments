function [X1, X2] = divides(X, alpha, verbose)
%========================================================================== 
% Syntax
%       [X1, X2] = divides(X, alpha, verbose);
%
% DIVIDES - Handles the divide phase (Approximate KNN)
%       The implementation of the algorithm is described in: "Fast 
%       Approximate KNN Graph Construction for High Dimensional Data via 
%       Recursive Lanczos Bisection" by Jie Chen, Haw-ren Fang, and Yousef
%       Saad.
%==========================================================================
% Inputs: 
%       X   - data vectors
%       alpha, verbose - (options ?)
%
% Outputs:   
%       X1, X2
%
% Uses L2distances
%==========================================================================
% Reference : "Fast Approximate KNN Graph Construction for High Dimensional
%              Data via Recursive Lanczos Bisection"
% Author   	: Avner Halevy
% Created	: Unknown
% Revised	: Dec 23, 2014 at 09:33 by Karamatou Yacoubou Djima
%==========================================================================

[~,n] = size(X);
if verbose; 
    fprintf('Divide! : |X|=%d',n); 
end

%Y = X([2,82],:); % least correlated bands
%Y = X([2,3],:); %2 arbitrary bands
Y = X(2:end,:); % all bands

c = mean(Y,2);
Y = Y - c*ones(1,n);

[~, ~, V]=svds(Y,1);
v = V(:,1);

sav = sort(abs(v));
h = sav(ceil(alpha*n));

ind1 = (v>=-h);
ind2 = (v<h);

%check if hyperplane failed
if sum(ind1)==0 || sum(ind2)==0
   h=floor(n/2);
   ind1 = 1:ceil(h + alpha*n);
   ind2 = floor(h-alpha*n):n;
end
X1 = X(:,ind1);
X2 = X(:,ind2);

if verbose; 
    fprintf(', |X1|=%d, |X2|=%d\n',sum(ind1),sum(ind2)); 
end


