% Uncomment each section to visualize the respective laplacian eigenmap

% Map for rings

G = gsp_ring(1000);
W = full(G.W); 
[mappedX, mapping, lambda] = lapbasic(W, 3, 1, 'JDQR');
for i=0:0.1:1
    plot(mappedX(:,1),mappedX(:,3))
    hold on
end

% % Map for spiral 
% 
% G = gsp_spiral(100,3);
% W = full(G.W);
% [mappedX, mapping, lambda] = lapbasic(W, 3, 1, 'JDQR');
% for i=0:0.1:1
%     plot3(((lambda(1))^(-i))*mappedX(:,1),((lambda(2))^(-i))*mappedX(:,2),((lambda(3))^(-i))*mappedX(:,3))
%     hold on
% end

% % Map for sphere
% 
% G = gsp_sphere(100);
% W = full(G.W);
% [mappedX, mapping, lambda] = lapbasic(W, 3, 1, 'JDQR');
% for i=0:0.1:1
%     plot3(((lambda(1))^(-i))*mappedX(:,1),((lambda(2))^(-i))*mappedX(:,2),((lambda(3))^(-i))*mappedX(:,3))
%     hold on
% end

% Map for swiss roll

% G = gsp_swiss_roll(500);
% W = full(G.W);
% [mappedX, mapping, lambda] = lapbasic(W, 3, 1, 'JDQR');
% % for i=0:0.1:1
% %     plot(((lambda(1))^(-i))*mappedX(:,1),((lambda(2))^(-i))*mappedX(:,2))
% %     hold on
% % end
% plot3(((lambda(1))^(-i))*mappedX(:,1),((lambda(2))^(-i))*mappedX(:,2),((lambda(3))^(-i))*mappedX(:,3),'o')

% Map for stochastic block graphs 

% G = gsp_stochastic_block_graph(1024,10);
% W = full(G.W);
% [mappedX, mapping, lambda] = lapbasic(W, 3, 1, 'JDQR');
% for i=0:0.1:1
%     plot3(((lambda(1))^(-i))*mappedX(:,1),((lambda(2))^(-i))*mappedX(:,2),((lambda(3))^(-i))*mappedX(:,3),'o')
%     hold on
% end