% Make a ring graph 
% 
N = 256;
G = gsp_spiral(N);

sub = 1:N/2;
% 
% % % A graph from an image 
% % 
% % Im = imread('Rcirc.png');
% % Im = imresize(Im,[32 32]);
% % Im = Im(:,:,1);
% % Im = im2double(Im);
% % G = im2graph(Im);
% % sub = 1:512;
% 
NDiff_Maps = Neumann_DiffMaps(G,sub,3,0.1);

[S, S_deltaS, S_UdeltaS, deltaS] = subs(G,sub);

Diff_Maps = gsp_laplacian_eigenmaps(S,3);

figure; 
subplot(2,2,1);
gsp_plot_graph(G);
title('Original Graph','interpreter','latex','FontSize',16);

subplot(2,2,2);
gsp_plot_graph(S);
title('Subgraph','interpreter','latex','FontSize',16);

subplot(2,2,3);
gsp_plot_graph(S_deltaS);
title('Subgraph with Boundary','interpreter','latex','FontSize',16);

subplot(2,2,4);
plot3(NDiff_Maps(:,3),NDiff_Maps(:,2),NDiff_Maps(:,1),'ro','DisplayName','Neumann');
hold on;
plot3(Diff_Maps(:,3),Diff_Maps(:,2),Diff_Maps(:,1),'bo','DisplayName','Diffusion');
title('Diffusion Embeddings','interpreter','latex','FontSize',16);
legend;

 
% % Compute the k largest (combinatorial) Neumann + Dirichlet spectra..
% k=10;
% [N,D, B, deltaT_S, T_S, N_mat] = Neumann_Dirichlet(G,sub);
% [V_N, lambda_N] = eigs(N,k,eps());
% [V_D, lambda_D] = eigs(D,k,eps()); 

% %Plot the spectra!!!!!!
% 
% % % For neumann, extend the values on the boundary via the N matrix.. 
% for i=1:k
%     f = N_mat*V_N(:,i);
%     gsp_plot_signal(S_deltaS,f);
%     pause(0.5);
% end
% % For dirichlet, extend the boundary by setting it zero!! 
% % 
% % g = vertcat(V_D(:,2),zeros(length(deltaS),1));
% % gsp_plot_signal(S_deltaS,g);

% NDiff_Maps = Neumann_DiffMaps(G,sub,3,1);

%  % Experiments with spheres 
%  size = 512;
%  sphere_graph = gsp_sphere(size);
%  coordinates = sphere_graph.coords;
%  elevation = coordinates(:,3);
%  polarcap = coordinates(elevation > 1/2,:);
%  %plot3(polarcap(:,1),polarcap(:,2),polarcap(:,3),'ro');
%  cap = find(elevation > 1/2);
%  distances = gsp_distanz(coordinates',coordinates');
%  eps=0.5;
%  weightmatrix = exp(-(1/(2*(eps)^2))*(distances.^2)) - eye(size);
%  S = graph(weightmatrix,'upper');
%  S = graph2struct(S);
%  S.coords = coordinates; 
%  
%  % Run Neumann Diffusion on the polar cap
%  NDiff_Maps = Neumann_DiffMaps(S,cap',5,1);
%  
%  % Run Standard Diffusion on the polar cap 
%  [T, T_deltaT, T_UdeltaT, deltaT] = subs(S,cap');
%  Diff_maps = My_Eigenmaps(T,1,5);
%  
%  %plot both
%  
%  figure; 
%  subplot(2,2,1);
%  plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'bo','DisplayName','Sphere');
%  hold on;
%  plot3(coordinates(cap,1), coordinates(cap,2), coordinates(cap,3),'ro','DisplayName','Cap');
%  title("Sphere and Polar Cap",'interpreter','latex','FontSize',16);
%  legend;
%  
%  subplot(2,2,2);
%  plot3(Diff_Maps(:,3),Diff_Maps(:,2),Diff_Maps(:,1),'ro');
%  title('Diffusion Map','interpreter','latex','FontSize',16);
%  
%  subplot(2,2,3);
%  plot(NDiff_Maps(:,1),NDiff_Maps(:,2),'bo');
%  title('2-D Neumann Map','interpreter','latex','FontSize',16);
%  
%  subplot(2,2,4);
%  plot3(NDiff_Maps(:,1),NDiff_Maps(:,2),NDiff_Maps(:,3),'bo');
%  title('3-D Neumann Map','interpreter','latex','FontSize',16);