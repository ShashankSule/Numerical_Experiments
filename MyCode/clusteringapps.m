%% Set up data here

%% Polar cap 

size = 1024;
param.sampling = 'uniform';
sphere_graph = gsp_sphere(size, param);

% set up the sphere
coordinates = sphere_graph.coords;
colorscore = 6*(coordinates(:,3) + 1) + ...
             (atan(coordinates(:,2)./coordinates(:,1))) + pi;
[~, inds] = sort(colorscore);
coordinates = coordinates(inds,:); 
cmap = parula(size); 
scatter3(coordinates(:,1), coordinates(:,2), coordinates(:,3), 25, ...
         cmap, 'filled');

% set up the cap 

elevation = coordinates(:,3);
polarcap = coordinates(elevation > 1/2,:);
cap = find(elevation > 1/2);
capcolour = parula(length(cap));
scatter3(polarcap(:,1), polarcap(:,2), polarcap(:,3), 25, capcolour, ...
         'filled');
%% Set up affinity matrix here

G = gsp_nn_graph(polarcap); 
K = G.W; 

H = gsp_nn_graph(coordinates);
[N,~, ~, ~, Ppi,~] = Neumann_Dirichlet(H,cap'); % Get the combinatorial 
                                               % Laplacian
K_neu = Ppi - N; % Adjacency matrix 
%% Pass to the diffusion solvers 
configParams.t = 1;
configParams.normalization = 'lb';
[NdiffMaps, N_Lambda, ~, ~, ~, ~] = calcDiffusionMap(K_neu, configParams);
[diffMaps, D_Lambda, ~, ~, ~, ~] = calcDiffusionMap(K, configParams);
LEMaps = gsp_laplacian_eigenmaps(G,3);
%% Plot! 

scatter3(LEMaps(:,1), LEMaps(:,2), LEMaps(:,3), 25, capcolour, ...
         'filled'); 
     
%% Clusters 

% Set up clusters here
nclusters = 5; 
std = 0.02; 
data = [normrnd(0,std, [100,1])  normrnd(0,std, [100,1]); %red 
           normrnd(1,std, [100,1])  normrnd(0,std, [100,1]); %blue 
           normrnd(0,std, [100,1])  normrnd(1,std, [100,1]); %green
           normrnd(-1,std, [100,1])  normrnd(0,std, [100,1]); %violet
           normrnd(0,std, [100,1])  normrnd(-1,std, [100,1])]; %yellow 
%names = ["red"; "blue"; "green", "violet", "yellow"];
colors = [repmat([1 0 0], 100, 1); 
          repmat([0 0 1], 100, 1); 
          repmat([0 1 0], 100, 1);
          repmat([255 0 255],100,1); 
          repmat([255 255 0], 100, 1)];
scatter(data(101:end,1), data(101:end,2), 65, colors(101:end, :), 'filled', ... 
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.0)

%% Set up affinity matrix here


clusterinds = 101:500; 
param.k = 10; 
param.sigma = 1; 
H = gsp_nn_graph(data,param);
[N,~, ~, ~, Ppi,~] = Neumann_Dirichlet(H,clusterinds); % Get the combinatorial 
                                               % Laplacian
K_neu = Ppi - N; % Adjacency matrix 
G = gsp_nn_graph(data(clusterinds,:), param); 
%% Pass to the diffusion solvers 
configParams.t = 1;
configParams.normalization = 'markov';
configParams.maxInd = 6; 
[NdiffMaps, N_Lambda, ~, ~, ~, ~] = calcDiffusionMap(K_neu, configParams);
%[NdiffMaps, N_Lambda, ~, ~, ~] = Neumann_diffusion_maps(data, ...
%                                   clusterinds, true, 4, 1.0); 

% Note that here calcDiffusionMap really just calculates the eigenvectors
% of K_neu, so running it for K_neu is the same as calculating the Neumann
% map and it should ideally be the same in Neumann_diffusion_maps but they
% seem to give WILDLY different answers and I am wont to believe Gal Mishne
% much more than myself. 

NdiffMaps = NdiffMaps(1:2, :); 
LEMaps_full = gsp_laplacian_eigenmaps(H,3);
LEMaps_sub = gsp_laplacian_eigenmaps(G,3); 
%NdiffMaps = [cos(pi/2) sin(pi/2); -sin(pi/2) cos(pi/2)]*NdiffMaps;
%% Plot! 
figure 

subplot(2,2,1)
scatter(data(101:end,1), data(101:end,2), 35, colors(101:end, :), 'filled', ... 
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.0)
title("Original")
subplot(2,2,2)
scatter(NdiffMaps(1,:), NdiffMaps(2,:), 35, colors(101:500, :), 'filled', ... 
       'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.0)
title("Neumann Maps"); 
subplot(2,2,3)
scatter(LEMaps_full(:,1), LEMaps_full(:,2), 35, colors, 'filled', ... 
         'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.0)
title("Laplacian Eigenmaps (entire data)"); 
subplot(2,2,4)
scatter(LEMaps_sub(:,1), LEMaps_sub(:,2), 35, colors(101:500, :), 'filled', ... 
         'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.0)
title("Laplacian Eigenmaps (unlabelled set)"); 
