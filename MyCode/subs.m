function [H,I,J, deltaS] = subs(G,sub)
% This function gives three important subgraphs: 
% H -- Subgraph of G induced by sub
% I -- Subgraph of G induced by sub with boundary edges 
% J -- Subgraph of G induced by sub union boundary vertices 
% deltaS -- the boundary
% G -- the graph
% sub -- selection of the vertices 

% First some mopping up. If G isn't a struct, we convert it to a struct
% object. If it is then we keep it. 

if ~isstruct(G)
    G = graph2struct(G);
end
% A stores the adjacency (or weight) matrix

A = G.W;

% find boundary vertices 

deltaS = boundary(G,sub); 

% find induced subgraph 

H = gsp_subgraph(G,sub);

% graph with subgraph and boundary 

% AdjS = zeros(length(sub)+length(deltaS)); % Creating an adjacency matrix
% AdjS(1:length(sub), 1:length(sub)) = A(sub,sub); 
% AdjS(length(sub)+1:end,1:length(sub)) = A(deltaS, sub);
% AdjS(1:length(sub), length(sub)+1:end) = (A(deltaS, sub))' ;
% I = graph(AdjS);
% I = graph2struct(I);

I = gsp_subgraph(G,[sub';deltaS]);
I.W(end-length(deltaS)+1:end,end-length(deltaS)+1:end) = zeros(length(deltaS));

% % graph induced by subgraph union boundary 

J = subg(G,[sub';deltaS]);

% AdjS_deltaS = A([sub reshape(deltaS, [1,length(deltaS)])], [sub reshape(deltaS, [1,length(deltaS)])]);
% % Creating reordered adjacency with subgraph vertices first
% J = graph(AdjS_deltaS);
% J = graph2struct(J);

end


