function [G, Conn_Comp] = nngraph(X, k, gamma, normalize)
%========================================================================== 
% Syntax
%       [G, Conn_Comp] = nngraph(X, k, gamma, normalize);
%
% NNGRAPH - Nearest Neighbors GRAPH
%       This routine computes a graph G using a Nearest Neighbors algorithm.
%==========================================================================
% Inputs: 
%       k       - number of dimensions.
%       gamma   - gamma = 0: spectral only, gamma = inf: spatial only,
%                 o/w uses fusion metric 
%       normalize - options.
%
% Outputs:   
%       G       - graph.
%       Conn_Comp - ?
%==========================================================================
% Reference : N/A
% Author   	: Unknown
% Created	: Unknown
% Revised	: Dec 23, 2014 at 10: by Karamatou Yacoubou Djima
%==========================================================================

% Size of X
sz = size(X);

if ndims(X) ~= 2
    X = reshape(X,[sz(1)*sz(2) sz(3)]);
end

%if ~exist('k', 'var')
 %   k = 12;
%end

X = X';

% Start graph computation
G = annfusion(X, sz, k, .1, 7500, gamma, 0);  
size(G)
G = max(G, G');
G = G .^ 2;
% G_before=G;
% if nargin >3
%     if add_spatial
%         if nargin == 4
%             load GS
%         end
%         G=max(G,GS);
%     end
% end

% Normalize graph
if normalize
    G = G./ max(max(G));
end
% G_before = G_before./max(max(G_before));

% Only embed largest connected component of the neighborhood graph
blocks = components(G)';
count = zeros(1, max(blocks));

for i = 1:max(blocks)
    count(i) = length(find(blocks == i));
end
[~, block_no] = max(count);

Conn_Comp = find(blocks == block_no);
% G = G(conn_comp, conn_comp);
% Conn_Comp = 1:length(G);

