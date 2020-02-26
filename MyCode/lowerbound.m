function G = lowerbound(A, lambda)
% Computes G(lambda) from  Lemma 0.1 
N = max(size(A));
G = 0 ;
for i = 1:N
    if sum(A(i,:)) < lambda
        G = G + 1 ;
    end 
end

