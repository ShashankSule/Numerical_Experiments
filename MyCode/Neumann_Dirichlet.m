function [N,D, B, deltaT_S, T_S, N_mat] = Neumann_Dirichlet(G,sub)
% Computes the Neumann and dirichlet operators 
% operator on the subgraph induced by sub
% G is the graph. Could be in struct type.  
% sub is the selection of vertices

if ~isstruct(G)
    G = graph2struct(G);
end

[~, ~, S_UdeltaS, ~] = subs(G,sub);
L = S_UdeltaS.L; % note that L is stored as sparse double
D = L(1:length(sub), 1:length(sub)); % the dirichlet matrix 
B = -L(length(sub)+1:end,1:length(sub)); %the boundary map
%l = S_deltaS.L;
v = B*ones(length(sub),1);
deltaT_S = diag(v); 
%deltaT_S = l(length(sub)+1:end, length(sub)+1:end);
N = sparse(D - (B')*(diag(1./diag(deltaT_S)))*B);
diagonal = diag(L);
T_S = diag(diagonal(1:length(sub)));
N_mat = vertcat(sparse(eye(length(sub))), sparse((diag(1./diag(deltaT_S)))*B));

end





