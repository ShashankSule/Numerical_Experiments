N = 10;
A = RandomGraph(N);
G = graph(A);
plot(G, 'NodeColor', 'r', 'MarkerSize', 10, 'LineWidth', 2);


% Path graphs

% vec =  zeros(1,10);
% vec(2) = 1;
% Adj = toeplitz(vec);
% G = graph(Adj)
% plot(G)

% Cyclic graph
% S
% vec = zeros(1,10);
% vec(2) = 1; vec(length(vec)) = 1; 
% Adj = toeplitz(vec);
% G = graph(Adj);
% plot(G)

% Complete graph
%weights = [1 -2 3 -4 5 -6 7 -8 9];
% vec = ones(1,9);
% vec(1) = 0;
% Adj = toeplitz(vec);

% Emb = embedding(N);
% G = graph(A);
% L = full(laplacian(G));
% [V, D] = eig(L);
%plot(G, 'XData', Emb(1,:), 'YData', Emb(2,:), 'ZData', V(:,6));
%stem(1:N, V(:,9))

