%Exploring edge crossings vs lower bound

N = 20;
figure
title('Comparing $$\mathcal{E}$$ to G');
for k = 1:9
    A = RandomGraph(N);
    G = graph(A);
    L = full(laplacian(G));
    [V, D] = eig(L);
    %eigs stores the eigenvalues 
    eigs = diag(D);
    %C stores the edge crossing number
    C = zeros(1,N);
    for i = 1:N
        C(i) = crossings(A, V(:,i));
    end
    %Creating the G vector. For some reason Matlab isn't allowing me to broad
    %cast so I'll just do it via a for loop. Sigh. 
    G = zeros(1,N);
    for i = 1:N
        G(i) = lowerbound(A, eigs(i));
    end
    subplot(3,3,k)
    plot(eigs(2:N), C(2:N), 'r.');
    hold on
    plot(eigs(2:N), G(2:N)/2, 'b-');
end

