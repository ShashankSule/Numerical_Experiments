% Script for computing a graph and the respective boundary 

% Compute the graph here. Some sample computations are given below 

% % Erdos Renyi random graph 
% 
% e = 50;
% v = 100;
% G := Graph::createRandomGraph(v,e, undirected):
% 
% Random Graph with a bernoulli distributed adjacency matrix
p=0.05;
N = 100;
A = RandomGraph(N,p);
G = graph(A);
plot(G)

% % Path graphs
% 
% vec =  zeros(1,10);
% vec(2) = 1;
% Adj = toeplitz(vec);
% G = graph(Adj);
% plot(G)
% 
% % Cyclic graph
% 
% vec = zeros(1,10);
% vec(2) = 1; vec(length(vec)) = 1; 
% Adj = toeplitz(vec);
% G = graph(Adj);
% plot(G)
% 
% % Complete graph
% weights = [1 -2 3 -4 5 -6 7 -8 9];
% vec = ones(1,9);
% vec(1) = 0;
% Adj = toeplitz(vec);

% Define the graph subset here in row or column vector form 

sub = randperm(100,20);

% create the three important subgraphs 
[S, S_deltaS, S_UdeltaS, deltaS] = subs(G,sub);

S = graph(S.W);
S_deltaS = graph(S_deltaS.W);
S_UdeltaS = graph(S_UdeltaS.W);
%deltaS = graph(deltaS.W)  

% define new colormap
mymap = [1 0 0
         0 0 0
         0 0 1];
colormap(mymap);

% plotting subgraph with boundary

colours = zeros(1,length(sub)+length(deltaS));
colours(1:length(sub)) = 1;
colours(length(sub)+1:end) = 2;

% subroutine for assigning colours to edges 

Edges = table2array(S_deltaS.Edges); 
Edges = Edges(:,1:2); % Extract the set of edges denoted by ordered pairs
e_colours = zeros(1,length(Edges)); % Set of edges 
for i = 1:length(e_colours)
    if Edges(i,1) <= length(sub) && Edges(i,2) <= length(sub)
        e_colours(i) = 1; % case when the edge is in S
    else
        e_colours(i) = 2;
    end
end

% end of subroutine

p = plot(S_deltaS, 'MarkerSize', 10, 'LineWidth', 2);
legend('Subgraph', 'Boundary');
p.NodeCData = colours;
p.NodeLabel = [];
p.EdgeCData = e_colours;
p


% % plotting graph induced by subgraph union boundary 
% 
% 
% v_colours = zeros(1,length(sub)+length(deltaS));
% v_colours(1:length(sub)) = 2;
% v_colours(length(sub)+1:end) = 3;
% 
% % % subroutine for assigning colours to edges 
% 
% Edges = table2array(S_UdeltaS.Edges); 
% Edges = Edges(:,1:2); % Extract the set of edges denoted by ordered pairs
% e_colours = zeros(1,length(Edges)); % Set of edges 
% for i = 1:length(e_colours)
%     if Edges(i,1) <= length(sub) && Edges(i,2) <= length(sub)
%         e_colours(i) = 2; % case when the edge is in S
%     elseif Edges(i,1) > length(sub) && Edges(i,1) > length(sub)
%         e_colours(i) = 3; % case when the edge is not in S or the boundary
%     else
%         e_colours(i) = 1;
%         
%     end
% end
% 
% % % end of subroutine
% 
% subgraph_location = embedding(length(sub),[0 0]);
% boundary_location = embedding(length(deltaS),[4 0]);
% location = [subgraph_location boundary_location];
% p = plot(S_UdeltaS, 'MarkerSize', 10, 'LineWidth', 2);
% legend('Subgraph', 'Boundary');
% p.NodeCData = v_colours;
% p.NodeLabel = [];
% p.EdgeCData = e_colours;
% p.XData = location(1,:);
% p.YData = location(2,:);
% p





