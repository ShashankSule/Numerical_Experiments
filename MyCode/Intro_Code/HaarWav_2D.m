im = imread('Rcirc.png');
im = imresize(im,[32 32]);
im = im(:,:,1);

[a,h,v,d] = haart2(im,'integer');
figure;

for i=0:8
    row = floor(i/3) + 1;
    column = mod(i,3) + 1; 
    subplot(3,3,i+1);
    imrec = ihaart2(a,h,v,d,i,'integer');
    colormap parula
    imagesc(imrec);
    title(strcat('Level', " ", num2str(8-i + 1)));
    axis off;
end

% to extract a 2^N x 2^N sized image, just pick d(N:end) and run ihaart

