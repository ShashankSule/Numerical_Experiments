function [affinitymatrix] = makeaffinitymat(mat)

%**************************************************************************
%   Detailed explanation goes here
% 2:45
% mat is a m x n matrix, m features, n observations
%**************************************************************************

[m, n] = size(mat);
m = m - 10;
affinitymatrix = zeros(m-10,n);

for i = 1:m
   for j = 1:m
      xi = mat(i,:);
      xj = mat(j,:);
      fi = isnan(xi);
      fj = isnan(xj);
      xi(fi) = 0;
      xj(fj) = 0; 
      ang_ij = 100*dot(xi,xj)/(norm(xi)*norm(xj));
      nan_ij = length(find(or(fi,fj)==1));
      aff_ij = ang_ij*(1 - nan_ij/n);
      affinitymatrix(i,j) = aff_ij;
   end  
end

end
