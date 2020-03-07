function [N,D] = Neumann_Dirichlet(G,sub)
% Computes the Neumann and dirichlet operators 
% operator on the subgraph induced by sub
% G is the graph. make sure it's in graph type! 
% sub is the selection of vertices 
[S, S_deltaS, S_UdeltaS, deltaS] = subs(G,sub);
L = laplacian(S_UdeltaS); % note that L is stored as sparse double
D = L(1:length(sub), 1:length(sub)); % the dirichlet matrix 
B = -L(length(sub)+1:end,1:length(sub)); %the boundary map
l = laplacian(S_deltaS);
deltaT_S = l(length(sub)+1:end, length(sub)+1:end);
N = sparse(D - (B')*(diag(1./diag(deltaT_S)))*B)
end





