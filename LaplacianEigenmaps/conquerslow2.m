function [G, D] = conquerslow2(G1, G2, D1, D2, k, verbose)
%========================================================================== 
% Syntax
%       [G, D] = conquerslow2(G1, G2, D1, D2, k, verbose);
%
% CONQUERSLOW2 - Handles the conquer phase (Approximate KNN)
%       The implementation of the algorithm is described in: "Fast 
%       Approximate KNN Graph Construction for High Dimensional Data via 
%       Recursive Lanczos Bisection" by Jie Chen, Haw-ren Fang, and Yousef
%       Saad.
%==========================================================================
% Inputs: 
%       G1, G2 - graphs
%       D1, D2 - distances
%       k      - numbers of k nearest neighbors
%       verbose- (?)
%
% Outputs:   
%       G	- graph
%       D  	- distances
%==========================================================================
% Reference : "Fast Approximate KNN Graph Construction for High Dimensional
%              Data via Recursive Lanczos Bisection"
% Author   	: Avner Halevy
% Created	: Unknown
% Revised	: Dec 23, 2014 at 09:23 by Karamatou Yacoubou Djima
%==========================================================================

if verbose; 
    fprintf('Conquer!: |G1|=%d, |G2|=%d, |G|=',size(G1,2),size(G2,2)); 
end

% Check if G1 and G2 are ordered
[G1,D1] = enforce_order(G1,D1);
[G2,D2] = enforce_order(G2,D2);

% Check G1 and G2 are fully unique
if numel(unique(G1(1,:))) < size(G1,2)
    [G1,D1] = enforce_unique(G1,D1,k);
end
if numel(unique(G2(1,:))) < size(G2,2)
    [G2,D2] = enforce_unique(G2,D2,k);
end

p1 = G1(1,:);
p2 = G2(1,:);

ov = intersect(p1,p2);
ovsz = numel(ov);

ov1 = ismember(p1,ov);
ov2 = ismember(p2,ov);

oG1 = G1(:,ov1);
oD1 = D1(:,ov1);

if sum(ov1) ~= ovsz || sum(ov2) ~= ovsz
    a = 3;
end

G1(:,ov1) = [];
D1(:,ov1) = [];

oG2 = G2(:,ov2);
oD2 = D2(:,ov2);
G2(:,ov2) = []; 
D2(:,ov2) = [];

% Why for loop
for i = 1:numel(ov)
    l = ismember(oG2(:,i),oG1(:,i));
    oD2(l,i) = inf;
end
oD = [oD1;oD2];
oG = [oG1;oG2];

[oDa, ii] = sort(oD(2:end,:),1);
oD(2:k+1,:) = oDa(1:k,:); 
oD = oD(1:k+1,:);

shift = (0:ovsz-1)*(2*k+1);
shift = repmat(shift,2*(k+1)-1,1);
try ii = shift+ii;
catch err
    a = 2;
end
oG=oG(1:k+1,:);

G=[G1 oG G2];
D=[D1 oD D2];
[s, p] = sort([G1(1,:) oG(1,:) G2(1,:)]);

G = G(:,p);
D = D(:,p);

if verbose; 
    fprintf('%d\n',size(G,2)); 
end

end

function [G,D] = enforce_order(G,D)

[~,I]=sort(G(1,:));
G=G(:,I);
D=D(:,I);
end

% Really bad coding, need to fix at some point
% Fixes problems arrising from non-unique input data
function [G,D]=enforce_unique(G,D,k)
i=1;
while i <= (size(G,2)-1)
    
    if G(1,i) == G(1,i+1)
        starti=i;
        i=i+1;
        while (i<=(size(G,2)-1)) && (G(1,i) == G(1,i+1))
            i=i+1;
        end
        finishi=i;        
        tempD = D(2:end,starti:finishi);
        tempG = G(2:end,starti:finishi);
        tempD = tempD(:);
        tempG = tempG(:);
        
        G(:,starti+1:finishi)=[];
        D(:,starti+1:finishi)=[];
               
        Ind=tempG~=G(1,starti);
        tempD=tempD(Ind);
        tempG=tempG(Ind);
        
        [B,I,J]=unique(tempG);
        
        tempG=tempG(I);
        tempD=tempD(I);
        
        [Y,I]=sort(tempD);
        
        G(2:end,starti)=tempG(I(1:k));
        D(2:end,starti)=tempD(I(1:k));
        
    else
       i=i+1; 
    end
end
end


