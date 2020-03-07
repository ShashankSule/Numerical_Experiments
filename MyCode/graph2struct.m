function [G] = graph2struct(H)
% Takes a graph type and outputs a graph struct to use for 
% the gsp box
A = adjacency(H);
N = max(size(A));
G.N = N;
G.W = A; 
% Assign a default circular embedding 
G.coords=[(cos((0:N-1)*(2*pi)/N))',(sin((0:N-1)*(2*pi)/N))'];
G.plotting.limits=[-1,1,-1,1];
G = gsp_graph_default_parameters(G);
end 

