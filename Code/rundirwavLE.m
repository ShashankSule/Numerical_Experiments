function [MapsX, Xrec, coeffs] = rundirwavLE(X)

% Just test

% X = imread(filename);
method = 'Shearlets';
sigma = 0.5;
scales = 3;
LE_Opts = [1, 30, 20, sigma];

[MapsX, Xrec, coeffs] = dirwavLE(X, method, scales, LE_Opts);

% Signal end of code    
load gong.mat; sound(y)
disp('Done!');

end