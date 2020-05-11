function [NDiff_maps] = Neumann_DiffMaps(G,sub,dim,t)
% G--Graph in struct form
% dim--number of dimensions
% subs--subset of vertices indexing subgraph
% t--scale

[N,~, ~, ~, T_S,~] = Neumann_Dirichlet(G,sub); % Extract the Neumann and degree matrices
n = length(sub);
R = eye(n) - ((T_S)^(-1/2))*N*((T_S)^(-1/2));
[A, S] = eigs(R,dim+1,'largestabs');
B = ((T_S)^(-1/2))*A;
NDiff_maps = B*(S^t);
NDiff_maps = NDiff_maps(:,2:end);
end





