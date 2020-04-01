figure;
N = 16;
data_matrix = zeros(1024,3,N);
for i=1:N
    % Generate plot
    z = randi([1 1024],1,1);
    M = zeros(32,32);
    col = rem(z-1,32) + 1;
    row = ceil(z/32);
    M(row, col) = 1;
    subplot(4,4,i);
    imagesc(M);
    
    %Make the data
    M = M';
    data_matrix(:,1,i) = M(:);
    for m=1:1024
            row = ceil(m/32);
            column = rem(m-1,32)+1;
            data_matrix(m,2,i) = row/32;
            data_matrix(m,3,i) = column/32;
    end
end
eps = 5;
weight_matrix = zeros(N,N);
for i = 1:N
    for j = 1:N
        u = (norm(data_matrix(:,:,i) - data_matrix(:,:,j)))^2;
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

