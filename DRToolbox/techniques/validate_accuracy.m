function p=validate_accuracy(data,training,validate)

%numel(training(:,1))
%size(data,2)

%if (numel(training(:,1))<size(data,2)+10)
%class=classify(data,data(training(:,1),:),training(:,2),'diagquadratic');
%class=classify(data,data(training(:,1),:),training(:,2),'mahalanobis');

%else
class=classify(data,data(training(:,1),:),training(:,2));
%end
class_check=class(validate(:,1));
p=sum(class_check==validate(:,2))/numel(validate(:,1))*100;
 
end