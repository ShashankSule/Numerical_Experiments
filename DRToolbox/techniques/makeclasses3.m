function [training,validate]=makeclasses3(data,mpix,npix,valfrac,mintrain)

%valfrac is the fraction of training data to be used for validation
%mintrain is the smallest number of training elements to allow for a class

training=[];
validate=[];

unique_class=unique(data(:,3));
count=1;
for i=1:numel(unique_class)
    
    ind=find(data(:,3)==unique_class(i));
    
    
    if numel(ind)-ceil(numel(ind)*valfrac)>mintrain
            data(ind,3)=count;
            count=count+1;

        I=randperm(numel(ind));
        validate = [validate ; ...
            data(ind( I( 1:ceil(numel(ind)*valfrac))),:)];
        training = [training ; ...
            data(ind( I( ceil(numel(ind)*valfrac)+1:end)),:)];
    end    
end

end
