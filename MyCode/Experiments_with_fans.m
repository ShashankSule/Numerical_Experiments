A = nfan(3);
G = graph(A);
D = diag(A*ones(max(size(A)),1));
plot(G);
L = full(laplacian(G));
L_norm = D^(-1/2)*L*D^(-1/2);
[V,D] = eig(L_norm);
