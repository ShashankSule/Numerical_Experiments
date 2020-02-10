function cropandrotate2(filename, nn, option)

% n size of image. 512x512 or 256x256. Must be power of 2

nsize = num2str(nn);

prefix = strcat('/Users/kyacouboudjima/Dropbox (Amherst College)/Research/Computing/Data/Placentae/Cropped',nsize);

img = imread(filename);

[n,m,~] = size(img);

nofr = floor(n/nn);
nofc = floor(m/nn);
imgnum = 1;

alpha = [90,180,270];

cd 

for i = 1:nofr
    for j = 1:nofc
        subimg = img(1+(i-1)*nn:i*nn, 1+(j-1)*nn:nn*j,:);
        % figure;imshow(subimg);
        subimgname = strcat(prefix,filename,'_Part_',num2str(imgnum),'_Angle_','0','.png');
        
        imwrite(subimg,subimgname,'png');
        if option == 1
            for k = 1:3
                sumimgrot = imrotate(subimg,alpha(k));
                % figure;imshow(sumimgrot);
                subimgname = strcat(prefix,filename,'_Part_',num2str(imgnum),'_Angle_',num2str(alpha(k)),'.png');
                imwrite(sumimgrot,subimgname,'png');
            end
        else
        end
        imgnum = imgnum + 1;
    end
end

end
