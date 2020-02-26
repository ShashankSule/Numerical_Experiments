function subg(G,sub)
% Find the graph induced by a subset of vertices 
% Arguments:
% 1. G is the ambient graph 
% 2. sub is a subset of {1,...,n} that identifies the vertices 
A = adjacency(G);
V_S = A(sub, sub);
S = graph(S);
return S
end

