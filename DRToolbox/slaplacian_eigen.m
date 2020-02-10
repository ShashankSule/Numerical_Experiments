function [mappedX, mapping] = slaplacian_eigen(X, no_dims, k, sigma, xpixs, ypixs, eig_impl)
%LAPLACIAN_EIGEN Performs non-linear dimensionality reduction using Laplacian Eigenmaps
%
%   [mappedX, mapping] = laplacian_eigen(X, no_dims, k, sigma, eig_impl)
%
% Performs non-linear dimensionality reduction using Laplacian Eigenmaps.
% The data is in matrix X, in which the rows are the observations and the
% columns the dimensions. The variable dim indicates the preferred amount
% of dimensions to retain (default = 2). The variable k is the number of 
% neighbours in the graph (default = 12).
% The reduced data is returned in the matrix mappedX.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.7.2b.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, 2010
% University California, San Diego / Delft University of Technology

xpixs=100;
ypixs=100;


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
	if ~exist('sigma', 'var')
		sigma = 1;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end
    
    Z=X(:,4);
    X=X(:,1:4);
    
    % Construct neighborhood graph
    disp('Constructing neighborhood graph...');
    if size(X, 1) < 4000
        G = L2_distance(X', X');
        
        % Compute neighbourhood graph
        [~, ind] = sort(G); 
        for i=1:size(G, 1)
            G(i, ind((2 + k):end, i)) = 0; 
        end
        G = sparse(double(G));
        G = max(G, G');             % Make sure distance matrix is symmetric
    else
        [G,ni] = find_nn(X, k);
    end
    G = G .^ 2;
    G = G ./ max(max(G));
    
    % Only embed largest connected component of the neighborhood graph
    blocks = components(G)';
    count = zeros(1, max(blocks));
    for i=1:max(blocks)
        count(i) = length(find(blocks == i));
    end
    [count, block_no] = max(count);
    conn_comp = find(blocks == block_no)
    %G = G(conn_comp, conn_comp);
    
    % Compute weights (W = G)
    disp('Computing weight matrices...');
    
%     H=zeros((size(G,1).^2-size(G,1))/2,1);
%     current=1;
%     for i=1:size(G,1)
%         
%         for j=i+1:size(G,1)
%            line = bresenham(reshape(1:100,10,10),[i j],0);
%            H(current)=norm(Z(line(2:end))-Z(line(1:end-1)));
%            current=current+1;
%         end
%     end
%     sigma_1=1;
%     H=squareform(H);
    %H=H(conn_comp,conn_comp);
    
    H=zeros(xpixs*ypixs,xpixs*ypixs);
    mat=reshape(1:xpixs*ypixs,xpixs,ypixs);
    size(ni);
    for i=1:xpixs*ypixs
         if mod(i,100)==0
               fprintf('%d percent done\n',100*(i/(xpixs*ypixs)));
           end
       for j=1:k
          %if (G(i,j)~=0)
          [x0, y0]=ind2sub(size(mat),i);
          [x1, y1]=ind2sub(size(mat),ni(i,j));
              %line = bresenham(mat,[i ni(i,j)]);
              line = bresenham(mat,[x0 x1;y0 y1],0);
              H(i,ni(i,j))=norm(Z(line(2:end))-Z(line(1:end-1)));
          %end
       end
    end
    G=G(conn_comp,conn_comp);
    H=H(conn_comp,conn_comp);
    H(H~=0)
    sigma_1=1;
    
    
    % Compute Gaussian kernel (heat kernel-based weights)
    G(G ~= 0) = exp(sigma_1*H(G~=0)).*exp(-G(G ~= 0) / (2 * sigma ^ 2));
        
    % Construct diagonal weight matrix
    D = diag(sum(G, 2));
    
    % Compute Laplacian
    L = D - G;
    L(isnan(L)) = 0; D(isnan(D)) = 0;
	L(isinf(L)) = 0; D(isinf(D)) = 0;
    
    % Construct eigenmaps (solve Ly = lambda*Dy)
    disp('Constructing Eigenmaps...');
	tol = 0;
    if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, lambda] = jdqz(L, D, no_dims + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, lambda] = eigs(L, D, no_dims + 1, tol, options);			% only need bottom (no_dims + 1) eigenvectors
    end
    
    % Sort eigenvectors in ascending order
    lambda = diag(lambda);
    [lambda, ind] = sort(lambda, 'ascend');
    lambda = lambda(2:no_dims + 1);
    
    % Final embedding
	mappedX = mappedX(:,ind(2:no_dims + 1));

    % Store data for out-of-sample extension
    mapping.K = G;
    mapping.vec = mappedX;
    mapping.val = lambda;
    mapping.X = X(conn_comp,:);
    mapping.sigma = sigma;
    mapping.k = k;
    mapping.conn_comp = conn_comp;
    
    