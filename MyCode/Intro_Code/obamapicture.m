im = imread('obama.jpg');
im = imresize(im, [512 512]);
[a,h,v,d] = haart2(im,'integer');
D = d(5:end);
H = h(5:end);
V = v(5:end);
Imz = ihaart2(a,H,V,D,1,'integer');
Imz = double(Imz(:,:,1));
G = im2graph(Imz);
S = gsp_compute_fourier_basis(G);
U = full(S.U);
fiedler_vector = U(:,2);
M = median(fiedler_vector);
classifier = fiedler_vector > 0; 
gsp_plot_signal(G,classifier);
colormap flag