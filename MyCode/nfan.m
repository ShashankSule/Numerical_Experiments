function A = nfan(N)
%Return an N fan graph with 2n+1 vertices
A = zeros(2*N+1); 
A(1,:) = ones(1,2*N+1); %The first vertex is the center of the fan
A(1,1) = 0;
for i = 1:N
    %i loops through each blade
    A(2*i,1) = 1;
    A(2*i + 1, 1) = 1;
    A(2*i, 2*i + 1) = 1;
    A(2*i + 1, 2*i) = 1;
end

end
    
    