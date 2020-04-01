% Laplacian eigenmaps with digital weights 

% First we provide a graph G. Enter graph here
N = 512; % Number of vertices 
G = gsp_spiral(N,7); % Creating a graph in struct version 

Diff_maps = My_Eigenmaps(G);

% Plotting the jdim diffusion map
%plot(Diff_maps(:,1),Diff_maps(:,2),'ro');
plot3(Diff_maps(:,1),Diff_maps(:,2),Diff_maps(:,3),'ro');

