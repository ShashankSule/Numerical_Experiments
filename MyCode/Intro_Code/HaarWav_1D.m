load BabyECGData;
% figure;
% p1 = plot(times,HR,'-');
% xlabel('Hours');
% ylabel('Heart Rate');
% p1.Color(4) = 0.25;
% hold on;
[a,d] = haart(HR,'integer');
% HaarHR = ihaart(a,d,1,'integer');
% plot(times,HaarHR,'Linewidth',1)
% title('Haar Approximation of Heart Rate')
imz = zeros(10,2048);

for i = 1:10
    HaarHR = ihaart(a,d,i,'integer');
    imz(i,:) = HaarHR';
end
figure;
subplot(2,1,1);
set(gca,'TickLabelInterpreter','latex');
set(groot, 'DefaultLegendInterpreter','latex');
p1 = plot(times,HR,'-','DisplayName', 'Heart Rate');
xlabel('Hours', 'interpreter','latex', 'FontSize',16);
ylabel('Heart Rate','interpreter', 'latex','FontSize',16); 
p1.Color(4) = 0.2;
hold on
p2 = plot(times,imz(3,:),'r-','LineWidth',1,'DisplayName','Level 8 Haar approximation');
hold on 
p3 = plot(times,imz(7,:),'k--','LineWidth',2,'DisplayName','Level 4 Haar approximation');
legend('FontSize',11);

subplot(2,1,2);
colormap copper
set(gca,'TickLabelInterpreter','latex');
set(groot, 'DefaultLegendInterpreter','latex');
image(imz,'CDataMapping', 'scaled');
cbh = colorbar;
cbh.Ticks = [];
ylabel(cbh, 'Heart Rate','interpreter','latex','FontSize',16);
xticks([]);
%xlabel("Hours", 'interpreter','latex','FontSize', 16);
ylabel("Scale (j)", 'interpreter','latex','FontSize', 16);

%plotting the haar functions

H = generate_haar(512);
times = linspace(0,1,512);
figure; 
for i=1:8
    if i==1
        titlestring = strcat("$\varphi$");
        minlim = -1;
        maxlim = 1;
    else
         j = floor(log2(i-1));
         k = (i-1)-2^j;
         titlestring = strcat("$\psi_{",num2str(j),",",num2str(k),"}$");
         maxlim = max(H(i,:));
         minlim = min(H(i,:));
    end
    subplot(2,4,i);
    set(gca,'TickLabelInterpreter','latex');
    set(groot, 'DefaultLegendInterpreter','latex');
    plot(times,H(i,:),'LineWidth',2,'DisplayName','Level 8 Haar approximation');
    yticks([]);
    xticks([0 0.25 0.50 0.75 1]);
    xticklabels({'0', '$\frac{1}{4}$', '$\frac{1}{2}$', '$\frac{3}{4}$', '1'});
    set(gca,'FontSize',16);
    legend('FontSize',11);
    title(titlestring,'interpreter','latex','FontSize',20);
    
    
end


