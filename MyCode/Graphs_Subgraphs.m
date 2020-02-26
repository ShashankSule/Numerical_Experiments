% script for computing a graph and the respective boundary 

% Compute the graph here. Some sample computations are given below 
% 
% % Erdos Renyi random graph 
% 
% e = 50;
% v = 100;
% G := Graph::createRandomGraph(v,e, undirected):



% Random Graph with a bernoulli distributed adjacency matrix
p=0.1;
N = 100;
A = RandomGraph(N,p);
G = graph(A);
plot(G)

% Path graphs

% vec =  zeros(1,10);
% vec(2) = 1;
% Adj = toeplitz(vec);
% G = graph(Adj)
% plot(G)

% Cyclic graph
% 
% vec = zeros(1,10);
% vec(2) = 1; vec(length(vec)) = 1; 
% Adj = toeplitz(vec);
% G = graph(Adj);
% plot(G)

% Complete graph
%weights = [1 -2 3 -4 5 -6 7 -8 9];
% vec = ones(1,9);
% vec(1) = 0;
% Adj = toeplitz(vec);

% Define the graph subset here in row or column vector form 

sub = randperm(100,25);

% create the corresponding subgraph

S = subg(G,sub);

% find boundary vertices 

deltaS = boundary(G,sub); 

% % graph with subgraph and boundary 

% AdjS = zeros(length(sub)+length(deltaS)); % Creating an adjacency matrix
% AdjS(1:length(sub), 1:length(sub)) = A(sub, sub); 
% AdjS(length(sub)+1:end,1:length(sub)) = A(deltaS, sub);
% AdjS(1:length(sub), length(sub)+1:end) = (A(deltaS, sub))' ;
% S_deltaS = graph(AdjS); 
% colours = zeros(1,length(sub)+length(deltaS));
% colours(1:length(sub)) = 2;
% colours(length(sub)+1:end) = 1;
% p = plot(S_deltaS, 'MarkerSize', 10, 'LineWidth', 2);
% p.NodeCData = colours;
% p.NodeLabel = [];
% p


% graph induced by subgraph union boundary 

AdjS_deltaS = A([sub reshape(deltaS, [1,length(deltaS)])], [sub reshape(deltaS, [1,length(deltaS)])]);
% Creating reordered adjacency with subgraph vertices first
S_UdeltaS = graph(AdjS_deltaS);
% Creating subgraph with boundary
v_colours = zeros(1,length(sub)+length(deltaS));
v_colours(1:length(sub)) = 2;
v_colours(length(sub)+1:end) = 1;

% % subroutine for assigning colours to edges 

Edges = table2array(S_UdeltaS.Edges); 
Edges = Edges(:,1:2); % Extract the set of edges denoted by ordered pairs
e_colours = zeros(1,length(Edges)); % Set of edges 
for i = 1:length(e_colours)
    if Edges(i,1) <= length(sub) && Edges(i,2) <= length(sub)
        e_colours(i) = 2; % case when the edge is in S
    elseif Edges(i,1) > length(sub) && Edges(i,1) > length(sub)
        e_colours(i) = 1; % case when the edge is not in S or the boundary
    else
        e_colours(i) = 0;
        
    end
end

% % end of subroutine
subgraph_location = embedding(length(sub),[0 0]);
boundary_location = embedding(length(deltaS),[4 0]);
location = [subgraph_location boundary_location];
p = plot(S_UdeltaS, 'MarkerSize', 10, 'LineWidth', 2);
p.NodeCData = v_colours;
p.NodeLabel = [];
p.EdgeCData = e_colours;
p.XData = location(1,:);
p.YData = location(2,:);
p





