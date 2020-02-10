function [training,validate]=makeclasses(data,mpix,npix,valfrac,mintrain)

%valfrac is the fraction of training data to be used for validation
%mintrain is the smallest number of training elements to allow for a class

training=[];
validate=[];

unique_class=unique(data(:,3));
for i=1:numel(unique_class)
    
    ind=find(data(:,3)==unique_class(i));
    if numel(ind)-ceil(numel(ind)*valfrac)>mintrain
        I=randperm(numel(ind));
        validate = [validate ; [sub2ind([mpix npix],...
            data(ind( I( 1:ceil(numel(ind)*valfrac))),1),...
            data(ind( I( 1:ceil(numel(ind)*valfrac))),2)) ...
            data(ind( I( 1:ceil(numel(ind)*valfrac))),3)] ];
        training = [training ; [sub2ind([mpix npix],...
            data(ind( I( ceil(numel(ind)*valfrac)+1:end)),1), ...
            data(ind( I( ceil(numel(ind)*valfrac)+1:end)),2)),...
            data(ind( I( ceil(numel(ind)*valfrac)+1:end)),3)] ];
    end    
end

end
