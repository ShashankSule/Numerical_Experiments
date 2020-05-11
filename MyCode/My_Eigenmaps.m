function [Diff_maps] = My_Eigenmaps(G,t,dim)
% Spits out 5 dimensional laplacian eigenmaps given a graph G in struct
% form
N = max(size(G.W));
% compute the random walk matrix here
% t is the scaling factor in diffusion map
lap = gsp_create_laplacian(G, 'normalized'); %lap will be a struct
D = diag(full(lap.W)*(ones(N,1))); %Stores diagonal matrix
M = eye(N) - full(lap.L); %Regularized random walk 
[X, Lambda] = eigs(M,dim+1,'largestabs'); %Compute spectral decomposition up to 5 eigenvectors
Phi = D^(1/2)*X; %Phi matrix
Psi = D^(-1/2)*X; %Psi matrix; you want to extract its columns! 
Diff_maps = (Psi)*(Lambda^t); %Multiplying each column with the respective eigenvalue
Diff_maps = Diff_maps(:,2:end); %Dropping the first column as it's all a constant
% Each column of diff_maps contains a coordinate of the diffusion map 
% to plot it, plot Diff_maps(:,j) in the jth coordinate
end

