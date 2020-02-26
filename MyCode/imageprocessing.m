rgb = im2double(imread('snapshot.png'));
%imshow(rgb)
[m,n,p] = size(rgb);
freq = fft2(rgb);
%image(abs(freq(:,:,3)))
red = rgb(:,:,1);
green = rgb(:,:,2);
blue = rgb(:,:,3);
red_dft = freq(:,:,1);
green_dft = freq(:,:,2);
blue_dft = freq(:,:,3);
abs_red_dft = abs(freq(:,:,1));
abs_green_dft = abs(freq(:,:,2));
abs_blue_dft = abs(freq(:,:,3));
%imshow(cat(3, recon(:,:,1), recon(:,:,2),zeros(m,n)))
% 
% %Plotting the various channels 
% subplot(2,2,1)
% imshow(rgb);
% title("Actual Image");
% subplot(2,2,2)
% imshow(cat(3, red,zeros(m,n),zeros(m,n)));
% title("Red Channel");
% subplot(2,2,3)
% imshow(cat(3,zeros(m,n),green,zeros(m,n)));
% title("Green Channel");
% subplot(2,2,4)
% imshow(cat(3,zeros(m,n),zeros(m,n),blue));
% title("Blue Channel");

% 
% % Plotting the DFT
% figure();
% subplot(2,2,1)
% image(abs(freq(:,:,1))), colorbar;
% title("DFT of the Red Channel");
% subplot(2,2,2)
% image(abs(freq(:,:,2)));
% title("DFT of the Green Channel");
% subplot(2,2,3)
% image(abs(freq(:,:,3)));
% title("DFT of the Blue Channel");
% subplot(2,2,4)
% imshow(rgb);
% title("Actual Image");
% 
% % Highpass filter: select only the bottom 99.99 frequencies
figure();

subplot(2,2,1)
filterval1 = prctile(abs(green_dft(:)),99);
green_dft(abs(green_dft) > filterval1) = 0;
filterval2 = prctile(abs(blue_dft(:)),99);
blue_dft(abs(blue_dft) > filterval2) = 0;
filterval3 = prctile(abs(red_dft(:)),99);
red_dft(abs(red_dft) > filterval3) = 0;
recon = ifft2(cat(3, red_dft,green_dft,blue_dft));
imshow(recon)
title("Filter using the top 1% of the DFT");

subplot(2,2,2)
freq = fft2(rgb);
red_dft = freq(:,:,1);
green_dft = freq(:,:,2);
blue_dft = freq(:,:,3);
filterval1 = prctile(abs(green_dft(:)),99.5);
green_dft(abs(green_dft) > filterval1) = 0;
filterval2 = prctile(abs(blue_dft(:)),99.9);
blue_dft(abs(blue_dft) > filterval2) = 0;
filterval3 = prctile(abs(red_dft(:)),99.9);
red_dft(abs(red_dft) > filterval3) = 0;
recon = ifft2(cat(3, red_dft,green_dft,blue_dft));
imshow(recon)
title("Filter using the top 0.1% of the DFT");

subplot(2,2,3)
freq = fft2(rgb);
red_dft = freq(:,:,1);
green_dft = freq(:,:,2);
blue_dft = freq(:,:,3);
filterval1 = prctile(abs(green_dft(:)),99.95);
green_dft(abs(green_dft) > filterval1) = 0;
filterval2 = prctile(abs(blue_dft(:)),99.99);
blue_dft(abs(blue_dft) > filterval2) = 0;
filterval3 = prctile(abs(red_dft(:)),99.99);
red_dft(abs(red_dft) > filterval3) = 0;
recon = ifft2(cat(3, red_dft,green_dft,blue_dft));
imshow(recon)
title("Filter using the top 0.01% of the DFT");

subplot(2,2,4)
freq = fft2(rgb);
red_dft = freq(:,:,1);
green_dft = freq(:,:,2);
blue_dft = freq(:,:,3);
filterval1 = prctile(abs(green_dft(:)),99.995);
green_dft(abs(green_dft) > filterval1) = 0;
filterval2 = prctile(abs(blue_dft(:)),99.999);
blue_dft(abs(blue_dft) > filterval2) = 0;
filterval3 = prctile(abs(red_dft(:)),99.999);
red_dft(abs(red_dft) > filterval3) = 0;
recon = ifft2(cat(3, red_dft,green_dft,blue_dft));
imshow(recon)
title("Filter using the top 0.001% of the DFT");

% Correlation plot

% x = randi([2 m-1],1,1000); y = randi([2 n-1],1,1000);
% 
% pixels = zeros(1,1000);
% neighbours = zeros(1,1000);
% 
% for i = 1:1000
%     pixels(i) = red(x(i),y(i));
%     neighbours(i) = red(x(i)+1, y(i) + 1);
% end
% 
% scatter(pixels,neighbours);
% xlabel("Pixel value");
% ylabel("Neighbour value");

% % Graph associated to the image 
% rgb = rand(10,10);
% n = 10;
% m = 10;
% vecvals = rgb(:);
% A = zeros(m*n);
% 
% for i = 1:m*n
%     for j = i+1:m*n
%         if(abs((vecvals(i)-vecvals(j))/vecvals(i)) <= 0.1)
%             A(i,j) = 1;
%         end
%     end
% end
% 
% Embedding = zeros(m*n, 2);
% for i = 1:m*n
%     c = ceil(i/n);
%     r = mod(i,m);
%     if(r==0)
%         r=m;
%     end
%     Embedding(i,1) = r;
%     Embedding(i,2) = c;
% end
% 
% A = A' + A;  
% G = graph(A);
% subplot(1,2,1)
% plot(G, 'XData', Embedding(:,2), 'YData', Embedding(:,1), 'ZData',vecvals);
% title("Graph of the image at 10% proximity");
% subplot(1,2,2)
% imshow(cat(3, rgb,zeros(m,n),zeros(m,n)));
% title("The image");