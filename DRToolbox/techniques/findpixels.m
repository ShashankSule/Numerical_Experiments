function new_pixelsClass = findpixels(pixelsClass,window)

%Takes training data for whole data set and window size and returns the
%training data that falls within the pixels considered (those with a full
%window around them). Also modifies their coordinates to match up with our
%dimension reduced data

window_radius = (window-1)/2;

new_pixelsClass = pixelsClass(find(pixelsClass(:,1)>window_radius...
    & pixelsClass(:,1)<201-window_radius...
    & pixelsClass(:,2)>window_radius...
    & pixelsClass(:,2)<201-window_radius),:);
for i=1:2
    new_pixelsClass(:,i)=new_pixelsClass(:,i)-window_radius;
end

end

