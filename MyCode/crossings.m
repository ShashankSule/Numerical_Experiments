function E = crossings(A,f)
E = 0;

N = max(size(A));
    for i = 1:N
        for j =i:N
            if (f(i)*f(j) < 0) && (A(i,j) == 1)
               E = E+1;
            end
        end
    end

end


