function startupLE(BaseDirectory)

%========================================================================== 
% startup_LE(BaseDirectory) adds various directories to the current path. 
% The code will assume that the current directory is the base directory 
% unless another directory is specified
%========================================================================== 

fprintf('startup_LE.m: setting LE paths ... \n');

if nargin==0
    Prefix  = [pwd filesep];
else
    Prefix  = [BaseDirectory filesep];
end

appendpath(([Prefix]));
appendpath(([Prefix 'Code']));
appendpath(([Prefix 'Data']));
appendpath(([Prefix 'DirectionalWavelets']));
appendpath(([Prefix 'DRToolbox']));
appendpath(([Prefix 'LaplacianEigenmaps']));

return

function appendpath(string)

fprintf('\t%s\\ \n', string);
p = genpath(string);
addpath(p);

return;

