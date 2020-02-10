function results=runall

%dim to test with
%dim=[1:10 15:5:50];
dim=1:50;
%load all window sizes
c=cell(1,5); 
load('Y_window1_knn12.mat');
c{1} = Y;
load('Y_window3_knn12.mat');
c{2} = Y;
load('Y_window5_knn12.mat');
c{3} = Y;
load('Y_window7_knn12.mat');
c{4} = Y;
load('Y_window9_knn12.mat');
c{5} = Y;

%get pixelsClass.mat
load('pixelsClass.mat')
pixelsClassSm = findpixels(pixelsClass,2*numel(c)-1);

p=3/4;

%hold the results
results=zeros(numel(c),numel(dim));

%number of runs to make
runs=25;

for i=1:runs
fprintf('Run #%d\n',i);
%loop over all training/validate sets
[training, validate] = makeclasses3(pixelsClassSm,200, 200, p,1);

for j=1:numel(c)
%loop over all cells
t=[ [ sub2ind( [200 - 2*(j-1), 200 - 2*(j-1)],...
    training(:,1) - (j-1), training(:,2)- (j-1))] training(:,3)]; 
v=[ [ sub2ind([200 - 2*(j-1),200 - 2*(j-1)],...
    validate(:,1) - (j-1), validate(:,2)- (j-1))] validate(:,3)];

for k=1:numel(dim)
%loop over all dim choices
results(j,k)=results(j,k)+validate_accuracy(c{j}(:,1:dim(k)),t,v);

end

end

end

results=results/runs;
% 
surf(dim,1:2:(2*numel(c)-1),results)
axis([1 dim(end) 1 (2*numel(c)-1) min(min(results)) 100])
xlabel('Dimensions');
ylabel('Window Size');
zlabel('Accuracy [%]');

end