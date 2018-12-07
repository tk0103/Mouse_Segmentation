%In = M1E2; InGT = M1GT; diff1 = 0.25; diff2 = 0.1;
%In = M2E2; InGT = M2GT; diff1 = 0.06; diff2 = 0.1;
In = M3E2; InGT = M3GT; diff1 = 0.06; diff2 = 0.07;

[N,edges] = histcounts(In(InGT ==1),'BinWidth',0.001);
[~,I] = max(N); modeval(1) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==2),'BinWidth',0.001);
[~,I] = max(N); modeval(2) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==3),'BinWidth',0.001);
[~,I] = max(N); modeval(3) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==4),'BinWidth',0.001);
[~,I] = max(N); modeval(4) = (edges(I) + edges(I+1)) /2;

Intemp = zeros(siz); mask_blaE2 = zeros(siz); Intemp(InGT == 1) = In(InGT == 1);
mask_blaE2(Intemp > (modeval(1) - diff1) & Intemp < (modeval(1) + diff1) ) = 1; 

Intemp = zeros(siz); mask_LkidE2 = zeros(siz); Intemp(InGT == 2) = In(InGT == 2);
mask_LkidE2(Intemp > (modeval(2) - diff2) & Intemp < (modeval(2) + diff2) ) = 1; 

Intemp = zeros(siz); mask_RkidE2 = zeros(siz); Intemp(InGT == 3) = In(InGT == 3);
mask_RkidE2(Intemp > (modeval(3) - diff2) & Intemp < (modeval(3) + diff2) ) = 1; 

cutM3GT = zeros(siz);  cutM3GT(InGT == 4) = 4;
mask_blaE2 = logical(mask_blaE2); cutM3GT(mask_blaE2) = 1; 
mask_LkidE2 = logical(mask_LkidE2); cutM3GT(mask_LkidE2) = 2;
mask_RkidE2 = logical(mask_RkidE2);cutM3GT(mask_RkidE2) = 3; 
cutM3GT(and(InGT == 1,not(mask_blaE2))) = 5;
cutM3GT(and(InGT == 2,not(mask_LkidE2))) = 6;
cutM3GT(and(InGT == 3,not(mask_RkidE2))) = 7;
cutM3GT = uint8(cutM3GT);

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
scatter3(ptemp_k1,ptemp_k2,ptemp_k3,'.');
hold on
scatter3(ptemp_b1,ptemp_b2,ptemp_b3,'.');
xlim([0.1 0.35]); xlabel('26-36 keV')
ylim([0.1 0.35]); ylabel('36-52 keV')
zlim([0.1 0.35]); zlabel('52-79 keV')
axis tight equal
%legend('left-kidney','background')