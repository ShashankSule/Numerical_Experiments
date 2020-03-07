function [G] = im2graph(I)
%takes an image as a matrix
%outputs a graph as a struct with lattice coords 

[n,m] = size(I);
W = zeros(n*m);
sigma_p = max(max(I));
coordinates = zeros(n*m,2);
for i=1:n*m
    row = floor((i-1)/n) + 1; 
    column = mod((i-1),m) + 1;
    coordinates(i,2) = 1 - (row/n) + (1/n);
    coordinates(i,1) = (column/m) - (1/m);
end

for i = 1:n*m 
    for j = 1:n*m
        row1 = floor((i-1)/n) + 1; %row position of the ith pixel
        column1 = mod((i-1),m) + 1; %column position of the ith pixel
        row2 = floor((j-1)/n) + 1; %row position of the ith pixel
        column2 = mod((j-1),m) + 1; %column position of the ith pixel
        g_distance = norm(coordinates(i,:)-coordinates(j,:));
        p_distance = I(row1, column1) - I(row2,column2);
        W(i,j) = exp(-((g_distance)^2)/(2))*exp(-((p_distance)^2)/(2*sigma_p));
    end
end
G.W = sparse(W);
G.coords = coordinates;
G.plotting.limits=[0,1,0,1];
G = gsp_graph_default_parameters(G);