function y = kernel_function_fast(v, X, center, kernel, param1, param2, type)
%KERNEL_FUNCTION Computes sum of (K * X) where X is a possible eigenvector
%
%   y = kernel_function(v, X, center, kernel, param1, param2)
%
% The function computes the sum of the elements of (K * v), where v is a
% possible eigenvector of K. This function is used to enable the use of
% EIGS in Kernel PCA. The other parameters of the function are the dataset
% X, the name of the kernel function (default = 'gauss'), and its
% corresponding parameters in param1 and param2.
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


global column_sums total_sum verbose

if verbose; fprintf('size(v)=[%d,%d]\n',size(v,1),size(v,2)); end

if ~exist('center', 'var')
    center = 0;
end
if ~exist('type', 'var')
    type = 'Normal';
end
if ~strcmp(type, 'ColumnSums'), fprintf(''); end

% If no kernel function is specified
if nargin == 2 || strcmp(kernel, 'none')
    kernel = 'linear';
end

% Construct result vector
y = zeros(size(X,2), size(v, 2));
n = size(X, 2);
if verbose; fprintf('In kernel_function (%d):\n',n); end

switch kernel
    
    case 'linear'
        % Retrieve information for centering of K
        if center || strcmp(type, 'ColumnSums')
            column_sums = zeros(1, n);
            for i=1:n
                % Compute single row of the kernel matrix
                K = X(:,i)' * X;
                column_sums = column_sums + K;
            end
            % Compute centering constant over entire kernel
            total_sum = ((1 / n^2) * sum(column_sums));
        end
        
        if ~strcmp(type, 'ColumnSums')
            % Compute product K*v
            for i=1:n
                % Compute single row of the kernel matrix
                K = X(:,i)' * X;
                
                % Center row of the kernel matrix
                if center
                    K = K - ((1 / n) .* column_sums) - ((1 / n) .* column_sums(i)) + total_sum;
                end
                
                % Compute sum of products
                y(i) = K * v;
            end
        else
            % Return column sums
            y = column_sums;
        end
        
    case 'poly'
        
        % Initialize some variables
        if ~exist('param1', 'var'), param1 = 1; param2 = 3; end
        
        % Retrieve information for centering of K
        if center || strcmp(type, 'ColumnSums')
            column_sums = zeros(1, n);
            for i=1:n
                % Compute column sums of the kernel matrix
                K = X(:,i)' * X;
                K = (K + param1) .^ param2;
                column_sums = column_sums + K;
            end
            % Compute centering constant over entire kernel
            total_sum = ((1 / n^2) * sum(column_sums));
        end
        
        if ~strcmp(type, 'ColumnSums')
            % Compute product K*v
            for i=1:n
                % Compute row of the kernel matrix
                K = X(:,i)' * X;
                K = (K + param1) .^ param2;
                
                % Center row of the kernel matrix
                if center
                    K = K - ((1 / n) .* column_sums) - ((1 / n) .* column_sums(i)) + total_sum;
                end
                
                % Compute sum of products
                y(i) = K * v;
            end
        else
            % Return column sums
            y = column_sums;
        end
        
    case 'gauss'
        tic
        try
        
        % Initialize some variables
        if ~exist('param1', 'var'), param1 = 1; end
        gap=750;

        % Retrieve information for centering of K
        
        if center || strcmp(type, 'ColumnSums')
            if length(column_sums) == 0
                fprintf('Calculating Column Sums:\n');
                column_sums = zeros(1, n);
                for i=1:gap:n
                    K = L2_distance(X(:,i:min(i+gap-1,n)), X);
                    K = exp(-(K.^2 / (2 * param1.^2)));
                    column_sums = column_sums + sum(K,1);
                    if verbose; fprintf('%d, ' ,i); end
                end
                if verbose; fprintf('\n'); end;
            end
            if length(total_sum) == 0
                % Compute centering constant over entire kernel
                total_sum = ((1 / n^2) * sum(column_sums));
            end
        end
        
        if ~strcmp(type, 'ColumnSums')
            % Compute product K*v
            
            fprintf('Calculating K*v\n');
            for i=1:gap:n
                K = L2_distance(X(:,i:min(i+gap-1,n)), X);
                K = exp(-(K.^2 / (2 * param1.^2)));
                if center
                    K = bsxfun(@minus,bsxfun(@minus,K, ((1 / n) .* column_sums)), ((1 / n) .* column_sums(i:min(i+gap-1,n)))') + total_sum;
                end
                
                
                y(i:min(i+gap-1,n),:)=K*v;
                
                
                if verbose; fprintf('%d, ' ,i); end;
                
            end
            if verbose; fprintf('\n'); end
            
            
        else
            % Return column sums
            y = column_sums;
        end
        
        catch err
           disp(err.message)
        end
        toc
    otherwise
        error('Unknown kernel function.');
end
%y=y;
end