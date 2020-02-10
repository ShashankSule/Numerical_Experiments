function thresholding(MapsX)

% Thresholding
Map = MapsX;
midthres = (max(max(Map)) + min(min(Map)))/2;
% midthres = 3*10^(-3);
nmax = 5;
nmin = 1;
jump = 2;
thresholdsneg = (-nmax:jump:-nmin)*midthres;
thresholdspos = (nmin:jump:nmax)*midthres;
thresholds = [thresholdsneg,thresholdspos];
%thresholds = thresholdsneg;
% thresholds = thresholdspos;

figure;
for j = 1:length(thresholds)
    Ind1 = Map > thresholds(j);
    Map(Ind1) = 1;
    subplot(2,length(thresholds)/2,j); colormap gray; imagesc(Map); colorbar;
    title(strcat('The Threshold is',{' '},num2str(thresholds(j))));
    Map = MapsX(:,:,1);
end

figure;
for j = 1:length(thresholds)
    Ind1 = Map > thresholds(j);
    Map(Ind1) = 1;
    subplot(2,length(thresholds)/2,j); colormap gray; imagesc(-Map); colorbar;
    title(strcat('The Threshold is',{' '},num2str(thresholds(j))));
    Map = MapsX(:,:,1);
end


figure;
for j = 1:length(thresholds)
    Ind1 = Map < thresholds(j);
    Map(Ind1) = 0;
    subplot(2,length(thresholds)/2,j); colormap gray; imagesc(Map); colorbar;
    title(strcat('The Threshold is',{' '},num2str(thresholds(j))));
    Map = MapsX(:,:,1);
end

figure;
for j = 1:length(thresholds)
    Ind1 = Map < thresholds(j);
    Map(Ind1) = 0;
    subplot(2,length(thresholds)/2,j); colormap gray; imagesc(-Map); colorbar;
    title(strcat('The Threshold is',{' '},num2str(thresholds(j))));
    Map = MapsX(:,:,1);
end