Img = zeros(32,32);
Img(16,16) = 1;
Img(16,17) = 1;
Img(17,16) = 1;
Img(17,17) = 1;


FT = fft2(Img);
[a,h,v,d] = haart2(Img, 'integer');

HT = cell2mat(d(1));

figure();


subplot(2,2,1);
imagesc(Img);
caxis('manual');
caxis([-1 1])
title('Original Image','interpreter','latex','FontSize',20);
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,2);
imagesc(HT);
caxis('manual');
caxis([-1 1])
title('16-point Haar Coefficients','interpreter','latex','FontSize',20);
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,3);
imagesc(real(FT));
caxis('manual');
caxis([-1 1])
title('Real part of Fourier coefficients','interpreter','latex','FontSize',20);
set(gca,'XColor', 'none','YColor','none')


subplot(2,2,4);
imagesc(imag(FT));
caxis('manual');
caxis([-1 1])
title('Imaginary part of Fourier coefficients','interpreter','latex','FontSize',20);
set(gca,'XColor', 'none','YColor','none');
cbh = colorbar;
cbh.Ticks = [-1 1];
ylabel(cbh, 'Luminescence/Coefficient Value','interpreter','latex','FontSize',20);


