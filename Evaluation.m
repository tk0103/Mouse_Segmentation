%%
%GraphCut Enrgy of each term
seikai = pM1GT(pmask1);
seikai = double(seikai);

Label = CurLabel;
term1 = RP((1:N)'+(Label(:)-1)*N);

Img = pM1E1(pmask1);
Sig = sub2ind(size(Sigmat), Label(GraphModel.Hi), Label(GraphModel.Hj));
Z = Img(GraphModel.Hi)-Img(GraphModel.Hj);
Zval = Z./Sigmat(Sig);
term2 = func(Zval,GraphModel.dist);
term2 = func(-Zval,GraphModel.dist);

shapeval = shape(Sig);
shapeval(shapeval == 2) = E1(shapeval == 2);
shapeval(shapeval == 3) = E2(shapeval == 3);
term3 = shapeval;

term1da = term1(Label == 3);
term2da = term2(Label == 3);
term3da = term3(Label == 3);

disp(sum(term1da(:))*lambda); 
disp(sum(term2da(:))./(1+c));
disp((sum(term3da(:))*c)./(1+c));
term = sum(term1da(:))+sum(term2da(:))+sum(term3da(:));
disp(term);

%%
Label = CurLabel;
E1test = pM1E1(pmask1);
P = GraphModel.Hi(Label(GraphModel.Hi)== 2 & Label(GraphModel.Hj)==4);
Q = GraphModel.Hj(Label(GraphModel.Hi)== 2 & Label(GraphModel.Hj)==4);

sig = 0.0174;
dif = (E1test(P) - E1test(Q))/(sig);
edges = [-2 -2:0.05:3 3];
xlim([-2 3]);
yyaxis left
histogram(dif,edges);
%histogram(dif);

x = edges';
z = exp(-x.^2 /2);
z(x < 0) = 1;

yyaxis right
xlabel('I_p - I_q / \sigma')
plot(x,z,'LineWidth',3);

%%
edge = 0:0.01:1.0;
hold on
histogram(pM3E2(pM3GT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(pM3E2(pM3GT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(pM3E2(pM3GT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(pM3E2(pM3GT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);
%xlabel('CT value')
%legend('Bladder','L','R','back')

%mutest =0.459;
%sigtest = sqrt(0.081);
mutest =0.4589;
sigtest = sqrt(0.0081);
y1 = pdf('Normal',edge,mutest,sigtest);
%y1 = y1./sum(y1(:));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

%mutest =0.350;
%sigtest = sqrt(0.0376);
mutest =0.35;
sigtest = sqrt(0.0376);
y2 = pdf('Normal',edge,mutest,sigtest);
%y2 = y2./sum(y2(:));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

%mutest =0.3775;
%sigtest = sqrt(0.0838);
mutest =0.3775;
sigtest = sqrt(0.0838);
y3 = pdf('Normal',edge,mutest,sigtest);
%y3 = y3./sum(y3(:));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

%mutest =0.2169;
%sigtest = sqrt(0.00074);
mutest =0.2169;
sigtest = sqrt(0.0007411);
y4 = pdf('Normal',edge,mutest,sigtest);
%y4 = y4./sum(y4(:));
plot(edge,y4,'Color',[204 153 255]/255,'LineWidth',2)
%%


%%
%energy2 
edge =[0 0.01:0.01:0.99 1.0];

bincounts = histc(pM3E2(pM3GT ==1),edge);
bin1 = bincounts./sum(bincounts(:))* (numel(pM3E2(pM3GT ==1))/ sumval);

bincounts = histc(pM3E2(pM3GT ==2),edge);
bin2 = bincounts./sum(bincounts(:))* (numel(pM3E2(pM3GT ==2))/ sumval);

bincounts = histc(pM3E2(pM3GT ==3),edge);
bin3 = bincounts./sum(bincounts(:))* (numel(pM3E2(pM3GT ==3))/ sumval);
sumval = numel(pM3E2(pM3GT ==1))+numel(pM3E2(pM3GT ==2))+numel(pM3E2(pM3GT ==3));


bar(edge,bin1,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[51 102 255]/255);
hold on
bar(edge,bin2,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 135 0]/255);
bar(edge,bin3,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 255 0]/255);


mutest =0.4589; sigtest = sqrt(0.0081);
y1 = pdf('Normal',edge,mutest,sigtest);
y1 = y1./sum(y1(:));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest =0.35; sigtest = sqrt(0.0376);
y2 = pdf('Normal',edge,mutest,sigtest);
y2 = y2./sum(y2(:));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest =0.3775; sigtest = sqrt(0.0838);
y3 = pdf('Normal',edge,mutest,sigtest);
y3 = y3./sum(y3(:));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)
xlim([0 1.0])
%%
%energy1 initial mouse2
edge =[0 0.01:0.01:0.99 1.0];
sumval = numel(pM1E2(pM1GT ==1))+numel(pM1E2(pM1GT ==2))+numel(pM1E2(pM1GT ==3));
bincounts = histc(pM1E2(pM1GT ==1),edge);
bin1 = bincounts./sum(bincounts(:))* (numel(pM1E2(pM1GT ==1))/ sumval);

bincounts = histc(pM1E2(pM1GT ==2),edge);
bin2 = bincounts./sum(bincounts(:))* (numel(pM1E2(pM1GT ==2))/ sumval);

bincounts = histc(pM1E2(pM1GT ==3),edge);
bin3 = bincounts./sum(bincounts(:))* (numel(pM1E2(pM1GT ==3))/ sumval);


bar(edge,bin1,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[51 102 255]/255);
hold on
bar(edge,bin2,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 135 0]/255);
bar(edge,bin3,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 255 0]/255);


mutest =1.047; sigtest = sqrt(0.78);
%mutest =0.487; sigtest = sqrt(0.0185);
y1 = pdf('Normal',edge,mutest,sigtest);
y1 = y1./sum(y1(:));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)


mutest =0.388; sigtest = sqrt(0.0425);
%mutest =0.334; sigtest = sqrt(0.0043);
y2 = pdf('Normal',edge,mutest,sigtest);
y2 = y2./sum(y2(:));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)


mutest =0.37; sigtest = sqrt(0.0188);
%mutest =0.352; sigtest = sqrt(0.0063);
y3 = pdf('Normal',edge,mutest,sigtest);
y3 = y3./sum(y3(:));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)
xlim([0 1.0])

%%
imagesc(pM1E2(:,:,206)');
axis tight equal off
colormap(gray)
caxis([ 0 0.7])
%%
%energy3
edge =[0 0.01:0.01:0.99 1.0];

bincounts = histc(pM3E3(pM3GT ==1),edge);
bin1 = bincounts./sum(bincounts(:))* (numel(pM3E3(pM3GT ==1))/ sumval);

bincounts = histc(pM3E3(pM3GT ==2),edge);
bin2 = bincounts./sum(bincounts(:))* (numel(pM3E3(pM3GT ==2))/ sumval);

bincounts = histc(pM3E3(pM3GT ==3),edge);
bin3 = bincounts./sum(bincounts(:))* (numel(pM3E3(pM3GT ==3))/ sumval);
sumval = numel(pM3E3(pM3GT ==1))+numel(pM3E3(pM3GT ==2))+numel(pM3E3(pM3GT ==3));

bar(edge,bin1,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[51 102 255]/255);
hold on
bar(edge,bin2,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 135 0]/255);
bar(edge,bin3,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 255 0]/255);

%mutest =0.3067; sigtest = sqrt(0.0022);
mutest =0.46; sigtest = sqrt(0.0787);
y1 = pdf('Normal',edge,mutest,sigtest);
y1 = y1./sum(y1(:));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

%mutest =0.3067; sigtest = sqrt(0.0022);
mutest =0.26; sigtest = sqrt(0.0058);
y2 = pdf('Normal',edge,mutest,sigtest);
y2 = y2./sum(y2(:));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

%mutest =0.268; sigtest = sqrt(0.0216);
mutest =0.265; sigtest = sqrt(0.003);
%y3 = pdf('Normal',edge,mutest,sigtest);
y3 = y3./sum(y3(:));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

%%
%energy4edge =[0 0.01:0.01:0.99 1.0];

bincounts = histc(pM3E4(pM3GT ==1),edge);
bin1 = bincounts./sum(bincounts(:))* (numel(pM3E4(pM3GT ==1))/ sumval);

bincounts = histc(pM3E4(pM3GT ==2),edge);
bin2 = bincounts./sum(bincounts(:))* (numel(pM3E4(pM3GT ==2))/ sumval);

bincounts = histc(pM3E4(pM3GT ==3),edge);
bin3 = bincounts./sum(bincounts(:))* (numel(pM3E4(pM3GT ==3))/ sumval);
sumval = numel(pM3E4(pM3GT ==1))+numel(pM3E4(pM3GT ==2))+numel(pM3E4(pM3GT ==3));

bar(edge,bin1,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[51 102 255]/255);
hold on
bar(edge,bin2,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 135 0]/255);
bar(edge,bin3,1.0,'EdgeColor',[50 50 50]/255,'FaceColor',[255 255 0]/255);

mutest =0.2549; sigtest = sqrt(0.0008);
y1 = pdf('Normal',edge,mutest,sigtest);
y1 = y1./sum(y1(:));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest =0.2344; sigtest = sqrt(0.0031);
y2 = pdf('Normal',edge,mutest,sigtest);
y2 = y2./sum(y2(:));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest =0.240; sigtest = sqrt(0.0066);
y3 = pdf('Normal',edge,mutest,sigtest);
y3 = y3./sum(y3(:));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

%%
%Distance feature
xim = zeros(siz2); yim = zeros(siz2);  zim = zeros(siz2); 
for i = 1:siz2(1)
    xim(i,:,:) = i;
end
for i = 1:siz2(2)
    yim(:,i,:) = i;
end
for i = 1:siz2(3)
    zim(:,:,i) = i;
end

%%
%Calculate kurtosis and skewness
class = 2;
kM1 = [kurtosis(pM1E1(pM1GT == class)); kurtosis(pM1E2(pM1GT == class)); kurtosis(pM1E3(pM1GT == class)); kurtosis(pM1E4(pM1GT == class));];
kM2 = [kurtosis(pM2E1(pM2GT == class)); kurtosis(pM2E2(pM2GT == class)); kurtosis(pM2E3(pM2GT == class)); kurtosis(pM2E4(pM2GT == class));];
kM3 = [kurtosis(pM3E1(pM3GT == class)); kurtosis(pM3E2(pM3GT == class)); kurtosis(pM3E3(pM3GT == class)); kurtosis(pM3E4(pM3GT == class));];

yM1 = [skewness(pM1E1(pM1GT == class)); skewness(pM1E2(pM1GT == class)); skewness(pM1E3(pM1GT == class)); skewness(pM1E4(pM1GT == class));];
yM2 = [skewness(pM2E1(pM2GT == class)); skewness(pM2E2(pM2GT == class)); skewness(pM2E3(pM2GT == class)); skewness(pM2E4(pM2GT == class));];
yM3 = [skewness(pM3E1(pM3GT == class)); skewness(pM3E2(pM3GT == class)); skewness(pM3E3(pM3GT == class)); skewness(pM3E4(pM3GT == class));];


%%
hold on
histogram(pM3E4(pM3GT==1),edge,'Facecolor','y','Normalization','probability','FaceAlpha',0.3);
%histogram(pM3E4(pM3GT==2),edge,'Facecolor','red','Normalization','probability','FaceAlpha',0.3);
%histogram(pM3E4(pM3GT==3),edge,'Facecolor','blue','Normalization','probability','FaceAlpha',0.3);
%%
hisval1 = pM2E1(pM2GT==2);  hisval2 = pM2E4(pM2GT==2);
hisval3 = pM2E1(pM2GT==2);  hisval4 = pM2E3(pM2GT==2);
hisval5 = pM2E1(pM2GT==2);  hisval6 = pM2E2(pM2GT==2);
hissize = 3000;
hisval1 = hisval1(1:hissize); hisval2 = hisval2(1:hissize);
hisval3 = hisval3(1:hissize); hisval4 = hisval4(1:hissize);
hisval5 = hisval5(1:hissize); hisval6 = hisval6(1:hissize);
%%
hold on
scatter(hisval1,hisval2,'.');
scatter(hisval3,hisval4,'.');
scatter(hisval5,hisval6,'.');
axis tight equal
xlim([0.15 0.4]);
ylim([0.15 0.4]);
xlabel('27-36 keV')
ylabel('36-52 keV, 52-79 keV, 79-118 keV')
%legend('E1&E4(¶t‘)','E1&E3(¶t‘)','E1&E2(¶t‘)')
%%
hold on
%scatter(pM3E1(pM3GT==3),pM3E4(pM3GT==3),'.');
%scatter(pM3E1(pM3GT==3),pM3E3(pM3GT==3),'.');
%scatter(pM3E1(pM3GT==3),pM3E2(pM3GT==3),'.');

scatter(pM2E1(pM2GT==1),pM2E4(pM2GT==1),'.');
scatter(pM2E1(pM2GT==1),pM2E3(pM2GT==1),'.');
scatter(pM2E1(pM2GT==1),pM2E2(pM2GT==1),'.');
axis tight equal
xlim([0.1 0.7]);
ylim([0.1 0.7]);
xlabel('27-36 keV')
ylabel('36-52 keV, 52-79 keV, 79-118 keV')
%legend('E1&E4(¶t‘)','E1&E3(¶t‘)','E1&E2(¶t‘)','E1&E4(δNγχ)','E1&E3(δNγχ)','E1&E2(δNγχ)','Location','northwest')
%%
class  = 1;
[h1,p1] = lillietest(pM1E3(pM1GT == class));
[h2,p2] = lillietest(pM2E3(pM2GT == class));
[h3,p3] = lillietest(pM3E3(pM3GT == class));

%%
%histogram
class = 2;
temp1 = M1E1(M1GT == class);
temp2 = M1E2(M1GT == class);
temp3 = M1E3(M1GT == class);
temp4 = M1E4(M1GT == class);


temp5 = wM1E1(M1GT == class);

%edges = [0 0:0.005:0.8 0.8];
edges = [0 0:0.01:1 1];
%edges = [-0.1 -0.1:0.005:0.5 0.5];

hold on
histogram(temp1,edges,'Normalization','probability');
histogram(temp2,edges,'Normalization','probability');
histogram(temp3,edges,'Normalization','probability');
histogram(temp4,edges,'Normalization','probability');
histogram(temp5,edges,'Normalization','probability');
%legend('26.9-36.0','36.0-51.9','51.9-78.9','78.9-119','Location','northeast')
%legend('Mouse1','Mouse2','Mouse3')
%%
st = 164; en = 439; 
pmask1 = zeros(siz); pmask2 = zeros(siz); pmask3 = zeros(siz);
pmask1(:,:,st:en) = mask1(:,:,st:en); pmask2(:,:,st:en) = mask2(:,:,st:en); pmask3(:,:,st:en) = mask3(:,:,st:en);
pmask1 = logical(pmask1); pmask2 = logical(pmask2); pmask3 = logical(pmask3); 
%%
%energyE1 = M1E1; energyE2 = M1E2; energyE3 = M1E3; energyE4 = M1E4; mask = pmask1;
%energyE1 = M2E1; energyE2 = M2E2; energyE3 = M2E3; energyE4 = M2E4; mask = pmask2;
energyE1 = M3E1; energyE2 = M3E2; energyE3 = M3E3; energyE4 = M3E4; mask = pmask3;
%%
IoutM3E1 = (energyE1 - prctile(energyE1(mask),3,1)) ./ (prctile(energyE1(mask),99.9,1) - prctile(energyE1(mask),3,1));
IoutM3E2 = (energyE2 - prctile(energyE2(mask),3,1)) ./ (prctile(energyE2(mask),99.9,1) - prctile(energyE2(mask),3,1));
IoutM3E3 = (energyE3 - prctile(energyE3(mask),3,1)) ./ (prctile(energyE3(mask),99.9,1) - prctile(energyE3(mask),3,1));
IoutM3E4 = (energyE4 - prctile(energyE4(mask),3,1)) ./ (prctile(energyE4(mask),99.9,1) - prctile(energyE4(mask),3,1));
%%
imagesc(M2E4(280:380,210:310,230)');
colormap(gray);
axis equal tight off
caxis([0 0.7]);
%%
temp = load_raw('C:\Users\yourb\Desktop\Animation\x64\Release\Data\M3GC2.raw','*uint8'); 
temp = reshape(temp,siz2);
temp(temp==4) = 0;
%%
imagesc(temp(255:365,235:345,66)');
colormap(map);
axis equal tight off
caxis([0 4]);
%%
%colormap gray
 rectangle('Position',[140,250,50,50],'FaceColor','none','EdgeColor','r',...
    'LineWidth',2)
%%
out = Imap;
out(Imap==4) =0;
%%

%%
slice = 66;
y =290;
rangex = [265,355];
rangey = [y-45,y+45];
test1 = pM3E2(:,y-15,slice);
test2 = pM3E2(:,y-15,slice);
%265:355,245:335

subplot(3,4,[1,6])
imagesc(pM3E2(:,:,slice)');
caxis([0,0.7])
colormap gray
axis equal tight off
rectangle('Position',[1,y-0.-15,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);


subplot(3,4,[3,8])
imagesc(pM3E2(:,:,slice)');
%colormap(map)
caxis([0,4])
axis equal tight off
 rectangle('Position',[1,y-0.5-15,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);

 
M1 = movmean(test1,3);
subplot(3,4,[9,10])
h= plot([test1]);
h(1).Color = 'r';
xlim(rangex);
ylim([0,0.7]);
%legend({'Gray value'});


M2 = movmean(test2,3);
subplot(3,4,[11,12])
h= plot([test2]);
h(1).Color = 'r';
xlim(rangex);
ylim([0,0.7]);
%legend({'Gray value'});

%%
imagesc(M1E2(130:220,240:330,360)');
colormap(gray);
axis equal tight off
caxis([0 0.7]);