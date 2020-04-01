figure;
N = 40;
data_matrix = zeros(256,N);
for i=1:N
    %z = randi([1 256],1,1);
    M = zeros(16,16);
    column = rem(i-1,16) + 1;
    M(7, column) = 1;
    data_matrix(:,i) = M(:);
    subplot(8,5,i);
    imagesc(M);
end
eps = 0.5;
weight_matrix = zeros(N,N);
for i = 1:N
    for j = 1:N
        u = (norm(data_matrix(:,i) - data_matrix(:,j)))^2;
        weight_matrix(i,j) = exp(-(1/(2*eps))*u);
    end
end

W = weight_matrix - eye(N);
G.N = N;
G.W = W;
% Assign a default circular embedding 
G.coords=[(cos((0:N-1)*(2*pi)/N))',(sin((0:N-1)*(2*pi)/N))'];
G.plotting.limits=[-1,1,-1,1];
G = gsp_graph_default_parameters(G);

