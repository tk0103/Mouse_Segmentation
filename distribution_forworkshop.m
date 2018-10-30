temp_k1 = pM3E1(pM3GT == 2);
temp_k2 = pM3E2(pM3GT == 2);
temp_k3 = pM3E3(pM3GT == 2);
temp_k4 = pM3E4(pM3GT == 2);
ptemp_k1 = temp_k1(1:1000);
ptemp_k2 = temp_k2(1:1000);
ptemp_k3 = temp_k3(1:1000);
ptemp_k4 = temp_k4(1:1000);
%%
temp_b1 = pM3E1(pM3GT == 4);
temp_b2 = pM3E2(pM3GT == 4);
temp_b3 = pM3E3(pM3GT == 4);
temp_b4 = pM3E4(pM3GT == 4);
ptemp_b1 = temp_b1(1:1000);
ptemp_b2 = temp_b2(1:1000);
ptemp_b3 = temp_b3(1:1000);
ptemp_b4 = temp_b4(1:1000);

%%
scatter3(ptemp_k1,ptemp_k2,ptemp_k3,'.');
hold on
scatter3(ptemp_b1,ptemp_b2,ptemp_b3,'.');
xlim([0.1 0.35]); xlabel('26-36 keV')
ylim([0.1 0.35]); ylabel('36-52 keV')
zlim([0.1 0.35]); zlabel('52-79 keV')
axis tight equal
%legend('left-kidney','background')
%%
edge = 0.01:0.0025:0.35;
hold on
histogram(pwM3E1(pM3GT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(pwM3E1(pM3GT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);






%%
%Mouse1, enegy2
edge =[0 0.01:0.01:2.99 3.0];

bincounts = histc(pM1E2(pM1GT ==1),edge);
bin1 = bincounts./sum(bincounts(:))* (numel(pM1E2(pM1GT ==1))/ sumval);

bincounts = histc(pM1E2(pM1GT ==2),edge);
bin2 = bincounts./sum(bincounts(:))* (numel(pM1E2(pM1GT ==2))/ sumval);

bincounts = histc(pM1E2(pM1GT ==3),edge);
bin3 = bincounts./sum(bincounts(:))* (numel(pM1E2(pM1GT ==3))/ sumval);
sumval = numel(pM1E2(pM1GT ==1))+numel(pM1E2(pM1GT ==2))+numel(pM1E2(pM1GT ==3));


bar(edge,bin1,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[51 102 255]/255);
hold on
bar(edge,bin2,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 135 0]/255);
bar(edge,bin3,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 255 0]/255);


mutest =0.772; sigtest = sqrt(0.311);
y1 = pdf('Normal',edge,mutest,sigtest);
y1 = y1./sum(y1(:));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest =0.362; sigtest = sqrt(0.0229);
y2 = pdf('Normal',edge,mutest,sigtest);
y2 = y2./sum(y2(:));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest =0.37; sigtest = sqrt(0.0118);
y3 = pdf('Normal',edge,mutest,sigtest);
y3 = y3./sum(y3(:));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)
xlim([0 3.0])