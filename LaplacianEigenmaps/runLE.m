function [mappedX, G] = runLE(IMG, ROI, k, no_dims, sigma, opt)

% options 2 x 1 vector:
% opt(1) = 0, keep entire image, else, use region of interest
% opt(2) = 0, LE for grayscale image, else, LE in each RGB band

% Retrieve data
% cd ~/Documents/MATLAB/RetDrusAnalysis
% [IMG, ~, ROI] = makepatientdata(patientname);

% Decide if you want smaller image
if opt(1) == 0 % if you keep entire image
    s = size(IMG{1});
    p1 = 1;
    p2 = s(1);
    q1 = 1;
    q2 = s(2);   
    m = s(1);
    n = s(2);

else
    % Use ROI (region of interest) to determine reduced image size
    p1 = ROI(1,1);
    p2 = ROI(1,2);
    q1 = ROI(2,1);
    q2 = ROI(2,2);
    m = length(p1:p2);
    n = length(q1:q2);
end

% Return in working directory
cd ~/Documents/MATLAB/LaplacianEigenmaps


if opt(2) == 0 % Compute embedding in grayscale
    % Initialize X
    X = zeros(m,n,9);
    for i = 1:9
        img = double(rgb2gray(IMG{i}(p1:p2,q1:q2,:)));
        X(:,:,i) = img;
    end
    
    % Make graph
    G = nngraph(X,k,0,1);
    
    % Compute embedding
    mappedX0 = lap_basic(G, no_dims, sigma, 'JDQR');   
    
    % Make final grayscale images
    for l = 1:no_dims
        mappedX(:,:,l) = reshape(mappedX0(:,l),[m,n]);
    end

else % Compute embedding in each color band RGB 
    % Initialize X
    X = zeros(m,n,3,9);

    % Reformat images for use in nngraph
    for i = 1:9
        for j = 1:3
            img = double(IMG{i}(p1:p2,q1:q2,j));
            X(:,:,i,j) = img;
        end
    end

    % Compute in each color band RGB
    % Make graph in each color band RGB
    G1 = nngraph(X(:,:,:,1),k,0,1);
    G2 = nngraph(X(:,:,:,2),k,0,1);
    G3 = nngraph(X(:,:,:,3),k,0,1);

    % Compute embedding in each color band RGB
    mappedX1 = lap_basic(G1, no_dims, sigma, 'JDQR');
    mappedX2 = lap_basic(G2, no_dims, sigma, 'JDQR');
    mappedX3 = lap_basic(G3, no_dims, sigma, 'JDQR');
    
    for l = 1:no_dims
        mappedX(:,:,1,l) = reshape(mappedX1(:,l),[m,n]);
        mappedX(:,:,2,l) = reshape(mappedX2(:,l),[m,n]);
        mappedX(:,:,3,l) = reshape(mappedX3(:,l),[m,n]);
    end
    
    G = [G1 G2 G3];
end


