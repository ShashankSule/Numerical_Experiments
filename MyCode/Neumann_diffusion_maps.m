function [NDmaps, Lambda, S, t, spectra] = Neumann_diffusion_maps(X, sub, autotune, dim, ...
   alpha, epsilon, delta, t)

% Compute neumann diffusions. Inputs: 
% X -- N x D data matrix 
% sub -- subset of data to be embedded
% autotune -- usually true
% dim -- dimension of the embedding
% alpha -- diffusion parameter 
% epsilon -- affinity tuning 
% delta -- scaling accuracy
% t -- scale 

% Sample test 

% %% Create data
% N = 128; % Number of vertices 
% param.sampling = 'uniform';
% G = gsp_sphere(N,param); % Creating a graph in struct version 
% 
% X = G.coords;
% sub = 1:N/2; 
% autotune = true;
% % delta = 0.01;
% alpha = 0.0; 
% dim = 3; 
% 
if autotune
    delta = 0.01;
end
%% First set up the affinity matrix from the data

n = size(X,1);
K = zeros(n);
Del = zeros(n);

%% Setting up the matrix of distances 

for i = 1:n
        for j = 1:n
                Del(i,j) = (norm(X(i,:)- X(j,:)))^2; 
        end
end

%% Set the epsilon parameter if on autotune
if autotune

    %Del = Del - diag(diag(Del)); 
    rowmins = zeros(n,1); 
    for i = 1:n 
        vec = Del(i,:);
        rowmins(i) = min(vec(vec > 0));
    end

    epsilon = 2*mean(rowmins);
end

%% Set up symmetric kernel
K = exp(-Del / (2*epsilon)); 

%% Normalize on alpha
q = K*ones(n,1); 
Q = diag(q); 
K_alph = (Q^(-alpha))*K*(Q^(-alpha));

%% Computing the Spectrum: Set up the matrices 
%K_alph = K_alph - diag(diag(K_alph)); 
H = graph(K_alph); % switch to the graph 
[N,~, ~, ~, Ppi,~] = Neumann_Dirichlet(H,sub); % extract neumann matrix
n = length(sub);
S = eye(n) - (Ppi^(-1/2))*N*(Ppi^(-1/2)); % Form random walk matrix
S = 0.5 * (S + S');
%% Computing the spectrum: Finally execute the algorithm
param.tol = 1e-6;
options.disp = 0;
options.isreal = 1;
options.issym = 1;
Afun = @(x) S*x;
S = sparse(S); 

[V, Lambda] = eigs(Afun, n, dim+1, 'lm', options); 

% [V,Lambda] = eigs(Afun, N, N, 'largestabs', ...
%                  'Tolerance', 1e-6,'IsFunctionSymmetric', true); 
% Compute the dim+1 largest eigenvalues and 
                           % eigenvectors
% [V, Lambda] = eig(S);       
[~,ind] = sort(diag(Lambda), 'descend'); % Sort the eigenvectors in ascending 
Lambda = Lambda(ind(1:dim+1), ind(1:dim+1));     % order 
V = V(:,ind(1:dim+1));

%% Computing the diffusion map 

% Normalizing once!
V        = V./ repmat(sqrt(sum(V.^2)),size(V,1),1);
Psi = sqrt(Ppi)*V; %The Psi matrix in the equation M = Phi*Lambda*Psi^T

%set scaling param

if autotune
    lambdas = diag(Lambda);
    t = ceil(log(delta)/(log(lambdas(end)) - log(lambdas(1))));
end

%t = 1;
% Set the eigenvalue tuning
Lambda_talpha = Lambda^(t-t*alpha);

% Compute the diff map

% Normalizing twice! 
NDmaps = (Psi./ repmat(sqrt(sum(Psi.^2)),size(Psi,1),1))*Lambda_talpha;

% Note here that the rows of Diff_maps are the points and columns the 
% spectrum. So we'll return the last dim columns 
spectra = eigs(S, size(S,1));
end