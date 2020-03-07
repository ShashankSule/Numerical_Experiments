% Make a ring graph

G = gsp_ring(100);

% Fix a subset of the vertices; say the first 50

sub = 1:50;

% create the important subgraphs 

[S, S_deltaS, S_UdeltaS, deltaS] = subs(G,sub); 

S = graph2struct(S); 

S_deltaS = graph2struct(S_deltaS); 

% Compute the k largest (combinatorial) Neumann + Dirichlet spectra..
k=6;
[N,D, N_mat] = Neumann_Dirichlet(G,sub);
[V_N, lambda_N] = eigs(N,k);
[V_D, lambda_D] = eigs(D,k); 

% Plot the spectra!!!!!!
% 
% % % For neumann, extend the values on the boundary via the N matrix.. 
for i=1:6
    f = N_mat*V_N(:,i);
    gsp_plot_signal(S_deltaS,f);
    pause(0.5);
end
% For dirichlet, extend the boundary by setting it zero!! 
% 
% g = vertcat(V_D(:,2),zeros(length(deltaS),1));
% gsp_plot_signal(S_deltaS,g);

