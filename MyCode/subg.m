function [S] = subg(G,sub)
% Find the graph induced by a subset of vertices 
% Arguments:
% 1. G is the ambient graph 
% 2. sub is a subset of {1,...,n} that identifies the vertices 
if ~isstruct(G)
    G = graph2struct(G);
end
S = gsp_subgraph(G,sub);
end

