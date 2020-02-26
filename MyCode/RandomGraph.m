function A = RandomGraph(N,p)

% Creates adjacency matrix for random graph with N vertices 
% Edge assigned between vertices i and j with probability p 
A = zeros(N,N);
for i = 1:N
    for j = i:N
        A(i,j) = binornd(1,p);
        A(j,i) = A(i,j);
    end
end

A = A - diag(diag(A));

end 
