% Laplacian eigenmaps with digital weights 

% First we provide a graph G. Enter graph here
N = 500; % Number of vertices 
G = gsp_spiral(N); % Creating a graph 

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
%plot(Diff_maps(:,1),Diff_maps(:,2),'ro');
plot3(Diff_maps(:,1),Diff_maps(:,2),Diff_maps(:,3),'ro');

% The neumann case 

% Define the graph subset here in row or column vector form 

sub = randperm(100,25);

% create the corresponding subgraph

S = subg(G,sub);

% find boundary vertices 

deltaS = boundary(G,sub); 

% % graph with subgraph and boundary 

AdjS = zeros(length(sub)+length(deltaS)); % Creating an adjacency matrix
AdjS(1:length(sub), 1:length(sub)) = A(sub, sub); 
AdjS(length(sub)+1:end,1:length(sub)) = A(deltaS, sub);
AdjS(1:length(sub), length(sub)+1:end) = (A(deltaS, sub))' ;
S_deltaS = graph(AdjS); 
% colours = zeros(1,length(sub)+length(deltaS));
% colours(1:length(sub)) = 2;
% colours(length(sub)+1:end) = 1;
% p = plot(S_deltaS, 'MarkerSize', 10, 'LineWidth', 2);
% p.NodeCData = colours;
% p.NodeLabel = [];
% p

% Creating the Neumann ma