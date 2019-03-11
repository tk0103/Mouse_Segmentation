%In = M1E2; InGT = M1GT; diff1 = 0.25; diff2 = 0.1;
In = M2E2; InGT = M2GT; diff1 = 0.06; diff2 = 0.1;
%In = M3E2; InGT = M3GT; diff1 = 0.06; diff2 = 0.07;
%In = M4E2; InGT = M4GT; diff1 = 0.06; diff2 = 0.16;

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

cutM2GT = zeros(siz);  cutM2GT(InGT == 4) = 4;
mask_blaE2 = logical(mask_blaE2); cutM2GT(mask_blaE2) = 1; 
mask_LkidE2 = logical(mask_LkidE2); cutM2GT(mask_LkidE2) = 2;
mask_RkidE2 = logical(mask_RkidE2);cutM2GT(mask_RkidE2) = 3; 
cutM2GT(and(InGT == 1,not(mask_blaE2))) = 5;
cutM2GT(and(InGT == 2,not(mask_LkidE2))) = 6;
cutM2GT(and(InGT == 3,not(mask_RkidE2))) = 7;
cutM2GT = uint8(cutM2GT);

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
Label = M;
E1test = M1E1(mask1);
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
slice = 73;
imagesc(M1E2(320:370,275:325,slice)');
axis tight equal off
caxis([0 4])
colormap(map)

hold on
temp = M1GT(320:370,275:325,slice);
temp(temp == 4) = 0;
v = [1,1];
contour(temp',v,'-r','LineWidth',2.0);
set(gca,'YDir','reverse')
%%
temp1 = M1E2(319:369,274:324,slice);
temp2 = M1E2(320:370,275:325,slice);
temp3 = abs(temp2 - temp1);
temp3 = exp(-temp3);
%%
imagesc(temp3');
colormap(gray)
axis tight equal off
caxis([0.8 1])

hold on
temp = M1GT(320:370,275:325,slice);
temp(temp == 4) = 0;
v = [1,1];
contour(temp',v,'-r','LineWidth',2.0);
set(gca,'YDir','reverse')
%%
slice = 70;
y =300;
rangex = [325,365];
rangey = [y-20,y+20];
test1 = M1E2(:,y,slice);
test2 = M1E2(:,y,slice);
%265:355,245:335

subplot(3,4,[1,6])
imagesc(M1GT(:,:,slice)');
caxis([0,4])
%colormap gray
colormap(map)
axis equal tight off
rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);

subplot(3,4,[3,8])
imagesc(M1E2(:,:,slice)');
%colormap(map)
%caxis([0,4])
axis equal tight off
 rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);
 
M1 = movmean(test1,3);
subplot(3,4,[9,10])
h= plot([test1]);
h(1).Color = 'r';
xlim(rangex);
ylim([0,3.0]);
%legend({'Gray value'});


M2 = movmean(test2,3);
subplot(3,4,[11,12])
h= plot([test2]);
h(1).Color = 'r';
xlim(rangex);
ylim([0,0.7]);
%legend({'Gray value'});

%%
GT = M3GT;
ptemp_k1 = M3E2(GT ==1);
ptemp_k2 = M3E3(GT ==1);
ptemp_k3 = M3E4(GT ==1);
ptemp_b1 = M3E2(GT ==4);
ptemp_b2 = M3E3(GT ==4);
ptemp_b3 = M3E4(GT ==4);

ptemp_k1 = ptemp_k1(1:3000);
ptemp_k2 = ptemp_k2(1:3000);
ptemp_k3 = ptemp_k3(1:3000);
ptemp_b1 = ptemp_b1(1:3000);
ptemp_b2 = ptemp_b2(1:3000);
ptemp_b3 = ptemp_b3(1:3000);
%%
scatter3(ptemp_k1,ptemp_k2,ptemp_k3,'.');
hold on
scatter3(ptemp_b1,ptemp_b2,ptemp_b3,'.');
xlim([0.1 0.5]); xlabel('26-36 keV')
ylim([0.1 0.5]); ylabel('36-52 keV')
zlim([0.1 0.5]); zlabel('52-79 keV')
axis tight equal
%legend('L.kidney','Background')
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 12;

%%
%edge =[0 0:0.025:2.8 2.8];
%xlim([0 2.8])
GT = M1GT;
edge =[0 0:0.005:0.6 0.6];
xlim([0 0.6])

hold on
histogram(wM3E1(M3GT == 1),edge,'Normalization','probability','EdgeAlpha',0.4);
histogram(wM3E1(M3GT == 4),edge,'Normalization','probability','EdgeAlpha',0.4);
%histogram(M3E2(M3GT == 1),edge,'Normalization','probability','EdgeAlpha',0.4);
%histogram(M4E2(M4GT == 1),edge,'Normalization','probability','EdgeAlpha',0.4);
legend('L.kidney','Background')
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 12;

%%
bn= [0.4867	,0.7922	,0.9534	,0.892]';
ln=	[0.6673	,0.7009	,0.7391	,0.7824]';
rn= [0.7921	,0.6464	,0.6308	,0.5112]';

bw=	[0.5467	,0.884	,0.8752	,0.367]';
lw= [0.5862	,0.5688	,0.4759	,0.3266]';
rw= [0.6832	,0.6918	,0.5163	,0.3114]';
%%
mb = [0.4867	,0.7922,	0.9534,	0.892]';
ml = [0.6673	,0.7009,	0.7391,	0.782]';
mr = [0.7921	,0.6464,	0.6308,	0.5112]';
mm = [0.6487	,0.713166667,	0.774433333,	0.7284]';

wb = [0.5467	,0.884	,0.8752	,0.367]';
wl = [0.5862	,0.5688	,0.4759	,0.3266]';
wr = [0.6832	,0.6918	,0.5163	,0.3114]';
wm = [0.605366667,	0.714866667	,0.622466667	,0.335]';

gb = [0.7471	,0.8777,	0.9258,	0.8304]';
gl = [0.79773161	,0.798019026	,0.218807725	,0.888751213]';
gr = [0.826829926	,0.75474252	,0.149646984	,0.783928965]';
gm = [0.790553845	,0.810153849	,0.431418236	,0.834360059]';

%%
boxplot([mb,wb,ml,wl,mr,wr,mm,wm],'Labels',{'Narrow(δNγχ)','Wide(δNγχ)'...
    ,'Narrow(¶t‘)','Wide(¶t‘)','Narrow(‰Et‘)','Wide(‰Et‘)','Narrow(•½‹Ο)','Wide(•½‹Ο)'});
ax = gca;
ax.FontSize = 10;
ylim([0 1.0])
ylabel('Jaccard Index')
%%
boxplot([mb,gb,ml,gl,mr,gr,mm,gm],'Labels',{'MAP(δNγχ)','GC(δNγχ)'...
    ,'MAP(¶t‘)','GC(¶t‘)','MAP(‰Et‘)','GC(‰Et‘)','MAP(•½‹Ο)','GC(•½‹Ο)'});
ax = gca;
ax.FontSize = 10;
ylim([0 1.0])
ylabel('Jaccard Index')
%%
[p1,h1] = signrank(bn,bw)

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
out = Imap2;
out(Imap2 == 4) = 0;
imagesc(out(100:290,150:340,200)');
%imagesc(out(:,:,200)');
axis tight equal off
caxis([0 4])
colormap(map);

%%
Mouse_voxnum = [38516 227013 216702 12316284; ...
                130557 220659 244271 14620045;...
                276949, 204858, 185717, 13149164;...
                248592 246458 235185 12864228];

c = categorical({'Mouse1','Mouse2','Mouse3','Mouse4'});
bar(c,Mouse_voxnum);
ax = gca;
ax.YScale = 'log';
ax.FontName = 'Times New Roman';
ax.FontSize = 15;
legend('Bladder','L.kidney','R.kidney','Background')
legend('boxoff')

%%
JI = [0.65615	0.7108	0.640525	0.68125	0.72795	0.780225	0.714675	0.792625	0.7515	0.651375	0.7551	0.754075	0.810675	0.8012	0.762075	0.679875;...
 0.535725	0.5876	0.4558	0.41425	0.69125	0.667725	0.53565	0.703525	0.60895	0.400225	0.717325	0.70005	0.663975	0.722325	0.725425	0.488975;...
 0.58005	0.613375	0.505425	0.413025	0.62505	0.67695	0.653375	0.666775	0.70385	0.496825	0.65215	0.623875	0.66675	0.645125	0.63295	0.547325;...
 0.590641667	0.637258333	0.533916667	0.502841667	0.681416667	0.7083	0.634566667	0.720975	0.6881	0.516141667	0.708191667	0.692666667	0.7138	0.722883333	0.706816667	0.572058333];

bar(JI',0.8,'EdgeColor','none');

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 12;
ax.XTickLabels = {'E1','E2','E3','E4','E1,2','E1,3','E1,4','E2,3','E2,4','E3,4','E1,2,3','E1,2,4','E1,3,4','E2,3,4','E1,2,3,4','Wide'};
ax.XTick = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; 
legend('Bladder','L.kidney','R.kidney','Background')
legend('boxoff')
ylim([0 1.0])

%%
imagesc(GT(:,:,slice)');
axis tight equal off
caxis([0 4])
colormap(map)

rectangle('Position',[110,150,180,180],'FaceColor','none','EdgeColor','r',...
    'LineWidth',2)
%270:380,210:320,
%%
slice = 215;
imagesc(GC(110:290,150:340,slice)');
axis tight equal off
caxis([0 4])
colormap(map)

hold on
temp = M3GT(110:290,150:340,slice);
temp(temp == 4) = 0;
v = [1,1];
contour(temp',v,'-r','LineWidth',2.0);
set(gca,'YDir','reverse')
axis tight equal off
%%
slice = 200;
imagesc(M1GT(110:300,155:345,slice)');
axis tight equal off
caxis([0 4])
colormap(map)

hold on
temp = GT(110:300,155:345,slice);
temp(temp == 4) = 0;
v = [1,1];
contour(temp',v,'-r','LineWidth',2.0);
set(gca,'YDir','reverse')

%%
Imap = load_raw('C:\Users\yourb\Desktop\new3\ImapM1.raw','*uint8');
%ImapW = load_raw('C:\Users\yourb\Desktop\new3\ImapWM4.raw','*uint8');
GC = load_raw('C:\Users\yourb\Desktop\new3\GCM1.raw','*uint8');
%GC2 = load_raw('C:\Users\yourb\Desktop\new3\GC2M4.raw','*uint8');
GT = M1GT;
Imap = reshape(Imap,siz); GC = reshape(GC,siz); 
%ImapW = reshape(ImapW,siz); GC2 = reshape(GC2,siz);
Imap(Imap == 4) = 0;  GT(GT == 4) = 0;  ImapW(ImapW == 4) = 0;
%%
JI = CalcuJI(GC2,M2GT,3);
disp(JI)


%%
In = M3E1(Rkidmask); InGT = M3GT(Rkidmask);
mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
%mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.01:0.7 0.7];
xlim([0 0.7])

hold on
%histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4,'FaceColor',[51 102 255]/255);
%histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4,'FaceColor',[255 135 0]/255);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);


mutest = mu(1,1); sigtest = sigma(1,1,1);
y1 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest = mu(2,1); sigtest =  sigma(1,1,2);
y2 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest = mu(3,1); sigtest = sigma(1,1,3);
y3 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)



%%
%325:365,280:320,70
Im = M4E2(150:230,230:310,198:202);
%Im = M3E2(140:220,255:335,198:202);
mask =  ones([81,81,5]); mask = logical(mask);
GraphModel2 = CreateFullyConnectedGraphWithMask(mask);
Z = (Im(GraphModel2.Hi)-Im(GraphModel2.Hj)).^2;

%%
val = GraphModel2.nodes;
%val = 1001;
n = 1;
out = zeros(val,1);
for n =1:val
   
    temp = Z(GraphModel2.Hi == n);
    temp2 = Z(GraphModel2.Hj == n);
    temp3 = GraphModel2.dist(GraphModel2.Hi == n);
    temp4 = GraphModel2.dist(GraphModel2.Hj == n);
    temp5 = vertcat(temp,temp2);
    temp6 = vertcat(temp3,temp4);
    val = exp(-temp5 ./(2*sigma(n)^2))./temp6;
    val2 = max(temp5);
    out(n,1) = min(val);
    out2(n,1) = min(val);
end

out = reshape(out,[81 81 5]);
%%
outexp = exp(-out ./2 ./0.0019);
%%
imagesc(out(:,:,3)');
axis tight equal off
colormap(gray)
caxis([0.58 0.60])

hold on
temp = M4GT(150:230,230:310,200);
%temp = M3GT(140:220,255:335,200);
temp(temp == 4) = 0;
v = [1,1];
contour(temp',v,'-r','LineWidth',2.0);
set(gca,'YDir','reverse')

%%
%325:365,280:320,70
Im = M1E2(325:365,280:320,68:72);
mask =  ones([41,41,5]); mask = logical(mask);
GraphModel2 = CreateFullyConnectedGraphWithMask(mask);
Z = (Im(GraphModel2.Hi)-Im(GraphModel2.Hj)).^2;

val = GraphModel2.nodes;
out = zeros(val,1);
out2 = zeros(val,1);
for n =1:val
    temp = Z(GraphModel2.Hi == n);
    temp2 = Z(GraphModel2.Hj == n);
    temp3 = GraphModel2.dist(GraphModel2.Hi == n);
    temp4 = GraphModel2.dist(GraphModel2.Hj == n);
    temp5 = vertcat(temp,temp2);
    temp6 = vertcat(temp3,temp4);
    
    val2 = max(temp5);
    val = exp(-val2 ./(2*0.05^2))./max(temp6);
    out(n,1) = val;
    out2(n,1) = val2;
end
out = reshape(out,[41 41 5]);
out2 = reshape(out2,[41 41 5]);
%%
imagesc(out(:,:,3)');
axis tight equal off
colormap(gray)
caxis([0.02 1.0])

hold on
temp = M1GT(325:365,280:320,70);
temp(temp == 4) = 0;
v = [1,1];
contour(temp',v,'-r','LineWidth',2.0);
set(gca,'YDir','reverse')

%%
[sigma,lambda] =ndgrid(0.0001:0.0004:0.008,0.001:0.002:0.07);
lambda = lambda(:); sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);

%%

[sigma,lambda] =ndgrid(0.0081:0.0004:0.01,0.001:0.002:0.07);
lambda = lambda(:); sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);

%%
[sigma,lambda,c] =ndgrid(0.0001:0.0004:0.008,0.001:0.002:0.07,0.94:0.02:0.98);
lambda = lambda(:); sigma = sigma(:); c = c(:);
OutputJI = zeros(size(sigma,1),1);

%%
[sigma,lambda,c] =ndgrid(0.0081:0.0004:0.01,0.001:0.002:0.07,0.94:0.02:0.98);
lambda = lambda(:); sigma = sigma(:); c = c(:);
OutputJI = zeros(size(sigma,1),1);
