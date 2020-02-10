function blocks = components(A)
%========================================================================== 
% Syntax
%       blocks = components(A);
%
% COMPONENTS - 
%       Finds connected components in a graph defined by the adjacency 
%       matrix A. The function outputs an n-vector of integers 1:k in 
%       blocks, meaning that A has k components. The vector blocks labels 
%       the vertices of A according to component. If the adjacency matrix A
%       is undirected (i.e. symmetric), the blocks are its connected compo-
%       nents. If the adjacency matrix A is directed (i.e. unsymmetric), 
%       the blocks are its strongly connected components.
%
% Notes: 
%       This file is part of the Matlab Toolbox for Dimensionality Reduct-
%       ion v0.7b. The toolbox can be obtained from 
%       http://ticc.uvt.nl/~lvdrmaaten
%       You are free to use, change, or redistribute this code in any way 
%       you want for non-commercial purposes. However, it is appreciated if
%       you maintain the name of the original author.
%
%==========================================================================
% Inputs: 
%       A   - adjacency matric
%
% Output:   
%       blocks - connected components in a graph defined by A
%
% See kknslow2
%==========================================================================
% Reference : N/A
% Author   	: (C) Laurens van der Maaten Tilburg University, 2008
% Created	: 2008
% Revised	: Dec 23, 2014 at 09:22 by Karamatou Yacoubou Djima
%==========================================================================

% Check size of adjacency matrix
[n, m] = size(A);
 if n ~= m
     error ('Adjacency matrix must be square')
 end
 
% Compute Dulmage-Mendelsohn permutation on A
if ~all(diag(A)) 
    [~, p, ~, r] = dmperm(A | speye(size(A)));
else
    [~, p, ~, r] = dmperm(A);  
end

% Compute sizes and number of clusters
sizes = diff(r);
k = length(sizes);

% Now compute the array blocks
blocks = zeros(1, n);
blocks(r(1:k)) = ones(1, k);
blocks = cumsum(blocks);

% Permute blocks so it maps vertices of A to components
blocks(p) = blocks;


