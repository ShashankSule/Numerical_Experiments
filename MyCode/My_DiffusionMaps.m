% Laplacian eigenmaps with digital weights 

% First we provide a graph G. Enter graph here
N = 128; % Number of vertices 
G = gsp_spiral(N,3); % Creating a graph in struct version 

Diff_maps = My_Eigenmaps(G,1);

% Plotting the jdim diffusion map
%plot(Diff_maps(:,1),Diff_maps(:,2),'ro');
plot3(Diff_maps(:,3),Diff_maps(:,4),Diff_maps(:,1),'ro');

