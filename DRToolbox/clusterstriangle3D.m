function [C] = clusterstriangle3D(n)

numClusters = 3;

% Define centers
c1 = [-2,0,0];
c2 = [2,0,0];
c3 = [0,0,12];

n = ceil(n/3);

shift1 = repmat(c1,[n,1]);
shift2 = repmat(c2,[n,1]);
shift3 = repmat(c3,[n,1]);


% Stretch the distance between the centers (optional)
%centers = 5*[c1; c2; c3];

% Compute the distance between centers
%D = L2_distance(centers', centers');
%minDistance = min(D(D > 0));
%n2 = n - (numClusters - 1) * 9;%
%k = 1;
%X = zeros([n,3]);

D = 3;
k = 3;

% Draw points on the unit sphere first (uniform angular)
X       = [randn(n,k),zeros(n,D-k)];
Xnorms  = sqrt(sum(X.^2,2));

for i=1:n
    X(i,:)=X(i,:)/Xnorms(i);
end

% Draw from data distribution for the radial density
beta=betarnd(k,1,[n 1]);

% Push the points in based on the radial distribution
for i=1:n
    X(i,:)=X(i,:)*beta(i);
end

% Clusters
X1 = X + shift1;
X2 = X + shift2;
X3 = X + shift3;

C = [X1;X2;X3];




            
            
            
            