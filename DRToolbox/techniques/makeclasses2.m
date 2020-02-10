function [training,validate]=makeclasses2(data,mpix,npix,numtrain)
%numtrain is the number of training elements to take for each class
%validate will be the remainder

training=[];
validate=[];

unique_class=unique(data(:,3));
for i=1:numel(unique_class)
    
    ind=find(data(:,3)==unique_class(i));
    if numel(ind)-2*numtrain>0
        I=randperm(numel(ind));
        training = [training ; [sub2ind([mpix npix],...
            data(ind( I( 1:numtrain)),1),...
            data(ind( I( 1:numtrain)),2)) ...
            data(ind( I( 1:numtrain)),3)] ];
        validate = [validate ; [sub2ind([mpix npix],...
            data(ind( I( numtrain+1:end)),1), ...
            data(ind( I( numtrain+1:end)),2)),...
            data(ind( I( numtrain+1:end)),3)] ];
    end    
end

end
