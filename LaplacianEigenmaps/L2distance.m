function D = L2distance(A, B)
%========================================================================== 
% Syntax
%       D = L2distance(A, B);
%
% L2DISTANCE - L2 (Euclidean distance)
%       This fully vectorized (VERY FAST!) routine computes the Euclidean 
%       distance between two vectors according to the equation
%                 ||A-B|| = sqrt (||A||^2 + ||B||^2 - 2*A.B )
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    D = L2distance(A,B);
%==========================================================================
% Inputs: 
%       A,B - (d x m) matrices between which distance is computed
%
% Outputs:   
%       D   - (n x n) Euclidean distances between vectors in A and B
%
%==========================================================================
% Reference : N/A
% Author   	: Roland Bunschoten, University of Amsterdam 
%             Intelligent Autonomous Systems (IAS) Group
%             Kruislaan 403  1098 SJ Amsterdam, bunschot@wins.uva.nl
%             (Tested on PC Matlab v5.2 and Solaris Matlab v5.3.)
% Revised   : JBT (3/18/00) to work for 1-dimensional vectors and to warn 
%             for imaginary numbers.  Also ensures that output is all real,
%             and allows the option of forcing diagonals to be zero. 
% Revised   : (C) Laurens van der Maaten Tilburg University, 2008
% Copyright : You are free to modify, extend and distribute this code 
%             granted that the author of the original code is mentioned as 
%             the original author of the code.
% Created	: Unknown
% Last Revised	: Dec 23, 2014 at 09:52 by Karamatou Yacoubou Djima
%==========================================================================

if nargin < 2
    error('Not enough input arguments');
end

if size(A, 1) ~= size(B, 1)
    error('A and B should be of same dimensionality');
end

if ~isreal(A) || ~isreal(B)
    warning('Computing distance table using imaginary inputs. Results may be off.'); 
end

% Padd zeros if necessray
if size(A, 1) == 1
    A = [A; zeros(1, size(A, 2))]; 
    B = [B; zeros(1, size(B, 2))]; 
end

% Compute distance table
D = sqrt(bsxfun(@plus, sum(A .* A)', bsxfun(@minus, sum(B .* B), 2 * A' * B)));

% Make sure result is real
D = real(D);

