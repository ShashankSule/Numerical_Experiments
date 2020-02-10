%========================================================================== 
% This script's goal is to play with the rungenLE code by changing the
% parameters, etc.
%========================================================================== 
% Written for Shashank's Thesis
%========================================================================== 

% First, set path using startup_LE. Use it in startup_LE form if your are 
% in the directory "NumericalExperiments." Otherwise, refer to the 
% startup_LE help file.)

% Declare variables and parameters (read their definition in rungenLE code)
I = imread('15082401 Placenta Fixed Barium Polarized.jpg');
k = 14; 
sigmaval = 1;
no_dims = 10;

% Format image for use in code
X = I(650:end-567,201:end-200,:); %This is quite large! Feel free to reduce
% the size as needed. 
X = X(1:256,1:256,:);
X = double(X);

% Call code
[MappedX, mapping, lambda] = rungenLE(X, k, no_dims, sigmaval); 