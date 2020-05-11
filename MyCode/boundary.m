function Bindex = boundary(G, sub)
% This function gives the boundary of a subgraph indexed by sub
% Arguments
% 1. G is the ambient graph encoded as a struct
% 2. sub indexes the subset 

if ~isstruct(G)
    G = graph2struct(G);
end

A = G.W; % Computing the adjacency matrix of G 
Bindex = zeros(length(sub),length(A)); % Stores the index values of the boundary vertices
for i=1:length(sub)
    for j=1:length(A(1,:))
        if A(sub(i),j) ~= 0 && isempty(sub(sub==j))
            Bindex(i,j) = j;
        end
    end
end

Bindex = unique(Bindex(:));
Bindex = Bindex(Bindex ~=0);

end 


    
    
    