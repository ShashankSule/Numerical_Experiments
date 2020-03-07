% Laplacian eigenmaps with digital weights 

% First we provide a graph G. Enter graph here
N = 512; % Number of vertices 
G = gsp_spiral(N,7); % Creating a graph 

% Now compute the random walk matrix here
t = 1; % scaling factor in diffusion map
lap = gsp_create_laplacian(G, 'normalized'); %lap will be a struct
D = diag(full(lap.W)*(ones(N,1))); %Stores diagonal matrix
M = eye(N) - full(lap.L); %Regularized random walk 
[X, Lambda] = jdqr(M); %Compute spectral decomposition up to 5 eigenvectors
Phi = D^(-1/2)*X; %This actually stores the diffusion maps in 5d
Diff_maps = (Phi)*Lambda; %Multiplying each column with the respective eigenvalue
Diff_maps = Diff_maps(:,2:end); %Dropping the first column as it's all a constant
% Each column of diff_maps contains a coordinate of the diffusion map 
% to plot it, plot Diff_maps(:,j) in the jth coordinate

% Plotting the jdim diffusion map
plot(Diff_maps(:,1),Diff_maps(:,2),'ro');
%plot3(Diff_maps(:,1),Diff_maps(:,2),Diff_maps(:,3),'ro');

