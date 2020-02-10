function [mappedX, mapping, lambda] = lapbasic(G, no_dims, sigma, eig_impl)
%========================================================================== 
% Syntax
%       [mappedX, mapping, lambda] = lapbasic(G, no_dims, sigma, eig_impl);
%
% LAPBASIC - Laplacian Basic
%       This routine computes the Laplacian Eigenmaps given a graph G.
%==========================================================================
% Inputs: 
%       G       - graph derived from a data matrix.
%       no_dims - number of dimensions.
%       sigma, eig_impl - options.
%
% Outputs:   
%       mappedX - Eigenmaps
%       mapping - (?)
%       lambda  - eigenvalues
%==========================================================================
% Reference : N/A
% Author   	: Unknown
% Created	: Unknown
% Revised	: Dec 23, 2014 at 10:17 by Karamatou Yacoubou Djima
%==========================================================================

% Checking nargin
if ~exist('no_dims', 'var')
    no_dims = 2;
end
if ~exist('sigma', 'var')
    sigma = 1;
end
if ~exist('eig_impl', 'var')
    eig_impl = 'Matlab';
end

% Compute weights (W = G)
% disp('Constructing L...');
 
% Compute Gaussian kernel (heat kernel-based weights)
G(G ~= 0) = exp((-G(G ~= 0) / (2 * sigma^2)));

% disp('Finished constructing Kernel');

% Construct diagonal weight matrix
D = diag(sum(G, 2));
% disp('Finished constructing D');

% Compute Laplacian
L = D - G;
% disp('Finished constructing L');

L(isnan(L)) = 0; D(isnan(D)) = 0;
L(isinf(L)) = 0; D(isinf(D)) = 0;

% Construct eigenmaps (solve Ly = lambda*Dy)
% disp('Constructing Eigenmaps...');
tol = 0;
if strcmp(eig_impl, 'JDQR')
    options.Disp = 0;
    options.LSolver = 'bicgstab';
    options.MaxIt=200; 
    % Only need bottom (no_dims + 1) eigenvectors
    %[MappedX, lambda] = eigs(diag(diag(D).^(-1/2))*L*diag(diag(D).^(-1/2)), no_dims + 1);	
    % Only need bottom (no_dims + 1) eigenvectors
    [mappedX, lambda] = jdqr(L, D, no_dims + 1, tol, options);			
    % Only need bottom (no_dims + 1) eigenvectors
    % [MappedXp, lambdap] = jdqz(Lp, D^-1, no_dims + 1, tol, options);			

elseif strcmp(eig_impl, 'SYM')
    options.Disp = 1;
    options.LSolver = 'bicgstab';
    options.MaxIt=200; 
    [mappedX, lambda] = eigs(diag(diag(D).^(-1/2))*L*diag(diag(D).^(-1/2)), no_dims + 1,'SM');
    % Only need bottom (no_dims + 1) eigenvectors 
    % [MappedX, lambda] = jdqr(diag(diag(D).^(-1/2))*L*diag(diag(D).^(-1/2)), no_dims + 1, tol, options);
    % Only need bottom (no_dims + 1) eigenvectors
    % [MappedXp, lambdap] = jdqz(Lp, D^-1, no_dims + 1, tol, options);			

elseif strcmp(eig_impl, 'UNNORM')
    options.Disp = 1;
    options.LSolver = 'bicgstab';
    options.MaxIt=200; 
    [mappedX, lambda] = eigs(L, no_dims + 1,'SM');
    % Only need bottom (no_dims + 1) eigenvectors
    % [MappedX, lambda] = jdqr(diag(diag(D).^(-1/2))*L*diag(diag(D).^(-1/2)), no_dims + 1, tol, options);			

elseif strcmp(eig_impl, 'eigifp')
    options.DISP=0;
    options.TOL=1e-7;
    [mappedX, lambda] = eigifp(L,D,no_dims+1,options);

elseif strcmp(eig_impl, 'jdcg')
    %options.jmax=50;
    options.Disp=1;
    options.v0=ones(size(G,1),1)+.1*rand(size(G,1),1);
    [mappedX, lambda]=jdcg_gep(L,D,no_dims+1,options);
    
elseif strcmp(eig_impl, 'schrod')
    alpha=.01;
    display(alpha);
    V = sparse(zeros(size(G,1),size(G,2)));
    V(430,430)=1;
    options.Disp = 1;
    options.LSolver = 'bicgstab';
    options.MaxIt=200; 
    [mappedX, lambda] = eigs(diag(diag(D).^(-1/2))*L*diag(diag(D).^(-1/2)) + alpha*V, no_dims + 1,'SM');
    % Only need bottom (no_dims + 1) eigenvectors
    %[MappedX, lambda] = jdqr(diag(diag(D).^(-1/2))*L*diag(diag(D).^(-1/2)), no_dims + 1, tol, options);	
    
else
    options.disp = 0;
    options.isreal = 1;
    options.issym = 1;
    [mappedX, lambda] = eigs(L, D, no_dims + 1, tol, options)	;		% only need bottom (no_dims + 1) eigenvectors
end

% Sort eigenvectors in ascending order
lambda = diag(lambda);
[lambda, ind] = sort(lambda, 'ascend');
no_dims=size(mappedX,2)-1;
lambda = lambda(2:no_dims + 1);

% Final embedding
mappedX = mappedX(:,ind(2:no_dims + 1));
% MappedX(:,1) = mappedX(:,1) / norm(mappedX(:,1));
% MappedX(:,2) = mappedX(:,2) / norm(mappedX(:,2));

% Store data for out-of-sample extension
mapping.K = G;
mapping.vec = mappedX;
mapping.val = lambda;
mapping.sigma = sigma;
