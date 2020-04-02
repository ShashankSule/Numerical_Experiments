% figure;
% N = 32;
% data_matrix = zeros(2,N);
% for i=1:N
%     % Generate plot
%     z = randi([1 32],1,1);
%     M = zeros(32,32);
%     col = rem(i-1,32) + 1;
%     
%     M(17, col) = 1;
%     subplot(4,8,i);
%     imagesc(M);
%     data_matrix(:,i) =  [17; col-16];
% %     %Make the data
% %     M = M';
% %     data_matrix(:,1,i) = M(:);
% %     for m=1:1024
% %             row = ceil(m/32);
% %             column = rem(m-1,32)+1;
% %             data_matrix(m,2,i) = row/32;
% %             data_matrix(m,3,i) = column/32;
%     end
% 
% 
% eps = 4;
% weight_matrix = zeros(N,N);
% for i = 1:N
%     for j = 1:N
%         u = (norm(abs(data_matrix(:,i)) - abs(data_matrix(:,j))));
%         weight_matrix(i,j) = exp(-(1/2)*(u/eps)^2);
%     end
% end
% 
% W = weight_matrix - eye(N);
% G.N = N;
% G.W = W;
% % Assign a default circular embedding 
% G.coords=[(cos((0:N-1)*(2*pi)/N))',(sin((0:N-1)*(2*pi)/N))'];
% G.plotting.limits=[-1,1,-1,1];
% G = gsp_graph_default_parameters(G);

% Embedding a torus 

% T = gsp_torus(32,32);
% coords = T.coords;
% Diff_maps = gsp_laplacian_eigenmaps(T,5);
% 
% figure;
% subplot(1,2,1);
% plot3(coords(:,1),coords(:,2),coords(:,3),'ro');
% title('Points sampled uniformly on a torus');
% 
% subplot(1,2,2);
% plot3(Diff_maps(:,2),Diff_maps(:,3),Diff_maps(:,1),'bo-');
% title('Diffusion Embedding of the torus as a pointcloud');

% Extracting the image 
Iz = imread('Rcirc.png');
im = imresize(Iz,[32 32]);
im = im(:,:,1);
H = im2graph(im2double(im));
t = [0 0.25 0.5 0.75 1];
for i=t
    Diff_maps = My_Eigenmaps(H,i);
    scale = num2str(i);
    labelstring = strcat('t=',scale);
    scatter3(Diff_maps(:,1),Diff_maps(:,2),Diff_maps(:,3),'DisplayName',labelstring);
    hold on;
end

legend;

%scatter(Diff_maps(:,1),Diff_maps(:,2),pointsize,Diff_maps(:,3))
%pointsize = 20; colorbar jet;
%scatter(Diff_maps(:,1),Diff_maps(:,2),pointsize,Diff_maps(:,3)); colormap jet;

% figure;
% subplot(2,2,1);
% imagesc(Iz(:,:,1));
% title('Original Image','interpreter','latex','FontSize',16);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% 
% subplot(2,2,2);
% imagesc(im);
% title('Pixellated Image', 'interpreter','latex','FontSize',16);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% 
% subplot(2,2,3);
% scatter3(Diff_maps(:,1),Diff_maps(:,2),Diff_maps(:,3));
% title('3-D Diffusion Embedding','interpreter','latex','FontSize',16);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% 
% subplot(2,2,4);
% pointsize = 20; 
% scatter(Diff_maps(:,1),Diff_maps(:,2),pointsize,Diff_maps(:,3)); colormap jet;
% title('Recreated Image','interpreter','latex','FontSize',16);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);