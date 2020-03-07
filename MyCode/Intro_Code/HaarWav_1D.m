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
p1 = plot(times,HR,'-','DisplayName', 'Heart Rate');
xlabel('Hours');
ylabel('Heart Rate'); 
p1.Color(4) = 0.2;
hold on
p2 = plot(times,imz(3,:),'-','LineWidth',2,'DisplayName','Level 8 Haar approximation');
hold on 
p3 = plot(times,imz(7,:),'-','LineWidth',2,'DisplayName','Level 4 Haar approximation');
legend()

subplot(2,1,2);
colormap copper
image(imz,'CDataMapping', 'scaled');
colorbar
xlabel("Time");
ylabel("Scale");

% %plotting the haar functions
% 
% H = generate_haar(512);
% times = linspace(0,1,512);
% figure; 
% for i=1:8
%     subplot(8,1,i);
%     plot(times,H(i,:));
%     yticks([-1 1]);
% end


