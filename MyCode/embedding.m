function E = embedding(N,k)
% N is the number of points
% k is the center of the circle 
n = 0:1:N-1;
x = cos(n*(2*pi/N)) + k(1);
y = sin(n*(2*pi/N)) + k(2);
E = [x;y]
end  

