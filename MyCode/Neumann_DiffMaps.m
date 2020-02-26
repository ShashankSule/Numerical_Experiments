% Enter the kind of graph here. Example is a ring

G = gsp_random_ring(100);

% Define the graph subset here in row or column vector form 

sub = randperm(100,25);

% create the corresponding subgraph

S = subg(G,sub);

% find boundary vertices 

deltaS = boundary(G,sub);

% graph with subgraph and boundary. code borrowed from graphs_subgraphs

AdjS = zeros(length(sub)+length(deltaS)); % Creating an adjacency matrix
AdjS(1:length(sub), 1:length(sub)) = A(sub, sub); 
AdjS(length(sub)+1:end,1:length(sub)) = A(deltaS, sub);
AdjS(1:length(sub), length(sub)+1:end) = (A(deltaS, sub))' ;
S_deltaS = graph(AdjS); 
% plot code
% colours = zeros(1,length(sub)+length(deltaS));
% colours(1:length(sub)) = 2;
% colours(length(sub)+1:end) = 1;
% p = plot(S_deltaS, 'MarkerSize', 10, 'LineWidth', 2);
% p.NodeCData = colours;
% p.NodeLabel = [];
% p


