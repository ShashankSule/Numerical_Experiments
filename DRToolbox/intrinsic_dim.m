function no_dims = intrinsic_dim(X, method,options)
%INTRINSIC_DIM Eestimate the intrinsic dimensionality of dataset X
%
%   no_dims = intrinsic_dim(X, method)
%
% Performs an estimation of the intrinsic dimensionality of dataset X based 
% on the method specified by method. Possible values for method are 'CorrDim'
% (based on correlation dimension), 'NearNbDim' (based on nearest neighbor 
% dimension), 'GMST' (based on the analysis of the geodesic minimum spanning
% tree), 'PackingNumbers' (based on the analysis of data packing numbers), 
% 'EigValue' (based on analysis of PCA eigenvalues), and 'MLE' (maximum 
% likelihood estimator). The default method is 'MLE'. All methods are
% parameterless.
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

  



    if ~exist('method', 'var')
        method = 'MLE';
    end
 
    if exist('options','var')
        if isfield('k','options'); k=options.k; else k=6; end
        if isfield('k1','options'); k1=options.k1; else k1=6; end
        if isfield('k2','options'); k2=options.k2; else k2=12; end
        if isfield('corrdimK','options'); corrdimK=options.corrdimK; else corrdimk=5; end

    else
        corrdimK=5;
        k=6;
        k1=6;
        k2=12;
    end
    
    % Handle PRTools dataset
    if strcmp(class(X), 'dataset')
        X = X.data;
    end
    
    % Remove duplicates from the dataset
    X = double(unique(X, 'rows'));

    % Make sure data is zero mean, unit variance
    X = X - repmat(mean(X, 1), [size(X, 1) 1]);
    X = X ./ repmat(var(X, 1) + 1e-7, [size(X, 1) 1]);
    
    switch method
        case 'CorrDim'
            % Compute correlation dimension estimation
            n = size(X, 1);
            %D = find_nn(X, 5);
            D=annfusion(X,corrdimK,.1,1000,0,false,false);
            
            
            %[foo, bar, val] = find(D);
            
            r1 = full(median(D(D>D(1,1)))); r2 = full(max(D(D>D(1,1))));
            s1 = 0; s2 = 0;
            X = X';
            XX = sum(X .^ 2);
            onez = ones(1,n);
            
            for i=1:n         
                p = X(:,i)';
                xx = sum(XX(:,i));
                xX = p * X;
                dist = xx * onez + XX - 2 * xX;
                dist = sqrt(dist(i+1:n));
                s1 = s1 + length(find(dist < r1));
                s2 = s2 + length(find(dist < r2));
                 if mod(n,1000) == 0
                         fprintf('%d, ',n)
                end
                
            end
            
            
            Cr1 = (2 / (n * (n - 1))) * s1;
            Cr2 = (2 / (n * (n - 1))) * s2;

            % Estimate intrinsic dimensionality
            no_dims = (log(Cr2) - log(Cr1)) / (log(r2) - log(r1));
            
        case 'NearNbDim'
            % Set neighborhood range to search in
            %k1 = 6;
            %k2 = 12;
            
            % Compute nearest neighbor dimension estimation
            [D,ind] = annfusion(X,k2,.1,1000,0,false,false);
            ind=ind(2:end,:)';
            Tk = zeros(1, k2 - k1);
            
            for k=k1:k2               
                fprintf('NearNbDim: %d of %d\n',k,k2)
                Tk(k - k1 + 1) = sum(full(D(sub2ind(size(D), (1:size(X, 1))', double(ind(:,k))))));
            end
            Tk = Tk ./ size(X, 1);
            no_dims = (log(Tk(end)) - log(Tk(1))) / (log(k2) - log(k1));
            
        case 'PackingNumbers'            
            % Parameters for the estimation
            r(1) = 0.1; r(2) = 0.5;
            epsilon = .01;
            max_iter = 20;
            done = 0;
            l = 0;
            
            % Perform iterations (until 'convergence')
            while ~done
                l = l + 1;
                perm = randperm(size(X, 1));
                
                % Compute L for two radiuses (size of C is packing number)
                for k=1:2
                   C = [];
                   for i=1:size(X, 1)
                       for j=1:numel(C)
                           if sqrt(sum((X(perm(i),:) - X(C(j),:)) .^ 2)) < r(k)
                               j = size(X, 1) + 1;
                               break;
                           end
                       end
                       if numel(C) == 0 || j < size(X, 1) + 1 
                           C = [C; perm(i)];
                       end
                   end
                   L(k, l) = log(numel(C));                 % maximum cardinality of an r(k)-separated subset of X
                end
                                
                % Estimate of intrinsic dimension
                no_dims = -((mean(L(2,:)) - mean(L(1,:))) / (log(r(2)) - log(r(1))));
                
                % Stop condition
                if l > 10
                    if 1.65 * (sqrt(var(L(1,:)) .^ 2 + var(L(2,:)) .^ 2) / (sqrt(l) * log(r(2)) - log(r(1)))) < no_dims * ((1 - epsilon) / 2)
                        done = 1;
                    end
                end
                if l > max_iter
                    done = 1;
                end
            end
            
        case 'GMST'         % Geodesic minimum spanning tree
            % Initialize some variables
            gamma = 1;
            M = 1; N = 10;
            samp_points = randperm(size(X,1));
            samp_points = samp_points(1:N);
            %k = 6;
            Q = length(samp_points);
            knnlenavg_vec = zeros(M, Q);
            knnlenstd_vec = zeros(M, Q);
            dvec = zeros(M, 1);
            
            % Compute Euclidean distance matrix
            %D = find_nn(X, k * 10); % wide range to deal with permutations
            D=annfusion(X,k*10,.1,1000,0,false,false);
            
            % Make M estimates
            for i=1:M
                
                % Perform resampling estimation of mean k-nn length
                j = 1;
                for n=samp_points
                    fprintf('Sample Point %d of %d\n',n,samp_points(end));
                    
                    % Sum cumulative distances over N random permutations
                    knnlen1 = 0;
                    knnlen2 = 0;
                    fprintf('Trial: ');
                    for trial=1:N
                        fprintf('%d ',trial);
                        
                        % Construct random permutation of data (throws out
                        % some points)
                        indices = randperm(size(X, 1));
                        indices = indices(1:n);
                        Dr = D(indices,:);
                        Drr = Dr(:,indices);
                        
                        % Compute sum of distances to k nearest neighbors
                        L = 0;
                        Drr = sort(Drr, 1);
                        for l=1:size(Drr, 2)
                            ind = min(find(Drr(:,l) ~= 0));
                            L = L + sum(Drr(ind + 1:min([ind + k size(Drr, 2)]), l));
                        end
                        
                        % Accumulate sum and squared sum over all trials
                        knnlen1 = knnlen1 + L;
                        knnlen2 = knnlen2 + L^2;
                    end
                    
                    % Compute average and standard deviation over N trials
                    knnlenavg_vec(i, j) = knnlen1 / N;
                    knnlenstd_vec(i, j) = sqrt((knnlen2 - (knnlen1 / N) ^ 2 * N) / (N - 1));

                    % Update counter
                    j = j + 1;
                    fprintf('\n');
                end

                % Compute least squares estimate of intrinsic dimensionality
                A = [log(samp_points)' ones(Q,1)];
                sol1 = inv(A' * A) * A' * log(knnlenavg_vec(i,:))';
                dvec(i) = gamma / (1 - sol1(1));
            end
            
            % Average over all M estimates
            no_dims = mean(abs(dvec));                
                        
        case 'EigValue'
            % Perform PCA
            [mappedX, mapping] = pca(X, size(X, 2));
            lambda = mapping.lambda ./ sum(mapping.lambda);
            
            % Plot eigenvalues
            plot(1:length(lambda), lambda)
            
            % Evaluate residuals
            no_dims = 0;
            while no_dims < size(X, 2) - 1 && lambda(no_dims + 1) > 0.025
                no_dims = no_dims + 1;
            end
            
        case 'MLE'
            % Set neighborhood range to search in
            %k1 = 6;
            %k2 = 12;

            % Compute matrix of log nearest neighbor distances
            %X = X';
            %[d n] = size(X);
            %X2 = sum(X.^2, 1); 
           % knnmatrix = zeros(k2, n);
            n=length(X);
           if n < 3000
                distance = repmat(X2, n, 1) + repmat(X2', 1, n) - 2 * X' * X;
                distance = sort(distance);
                knnmatrix= .5 * log(distance(2:k2 + 1,:));
            else
%                 for i=1:n
%                     distance = sort(repmat(X2(i), 1, n) + X2 - 2 * X(:,i)' * X);
%                     distance = sort(distance);
%                     knnmatrix(:,i) = .5 * log(distance(2:k2 + 1))'; 
%                     if mod(n,1000) == 0
%                         fprintf('%d, ',n)
%                     end
%                 end

                    
                D=annfusion(X,k2,.1,1000,0,false,false);
                D=sort(full(reshape(D(D>D(1,1)),k2,n)));
                knnmatrix = .5 * log(D); 

            end  
            % Compute the ML estimate
            S = cumsum(knnmatrix, 1);
            indexk = repmat((k1:k2)', 1, n);
            dhat = -(indexk - 2) ./ (S(k1:k2,:) - knnmatrix(k1:k2,:) .* indexk);

            % Plot histogram of estimates for all datapoints
            %hist(mean(dhat), 80), pause
            
            % Average over estimates and over values of k
            no_dims = mean(mean(dhat));
            
        otherwise
            error('Unknown method for estimating intrinsic dimensionalities.');
    end
    
