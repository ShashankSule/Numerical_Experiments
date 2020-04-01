% Let's create a circulant matrix! 

i = 10;
v = [0 ones(1,i) zeros(1,99-(2*i)) ones(1,i)];

A = toeplitz([v(1) fliplr(v(2:end))], v);
 

P = eye(100); 
P = P(randperm(100),:);

Adj = P*A*P';

G = graph(Adj);

G = graph2struct(G);

Diff_Maps = My_Eigenmaps(G);

% Plotting the jdim diffusion map
subplot(1,2,1)
colormap gray
imagesc(Adj);
title("Adjacency Matrix of a circulant graph $V=100$, $k=20$",'fontsize',16,'interpreter','latex');

subplot(1,2,2);
scatter(Diff_maps(:,1),Diff_maps(:,2));
title("2 Dimensional Diffusion embedding",'fontsize',16,'interpreter','latex');
%plot3(Diff_maps(:,1),Diff_maps(:,2),Diff_maps(:,3),'ro');




