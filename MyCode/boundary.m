function Bindex = boundary(G, sub)
% This function gives the boundary of a subgraph indexed by sub
% Arguments
% 1. G is the ambient graph
% 2. sub indexes the subset 

A = adjacency(G); % Computing the adjacency matrix of G 
Bindex = zeros(length(sub),length(A)); % Stores the index values of the boundary vertices
for i=1:length(sub)
    for j=1:length(A(1,:))
        if A(sub(i),j) == 1 && isempty(sub(sub==j))
            Bindex(i,j) = j;
        end
    end
end

Bindex = unique(Bindex(:));
Bindex = Bindex(Bindex ~=0);

end 


    
    
    