bn= [0.5705,0.794,0.9536,0.8956]';
ln=	[0.6686,0.7118,0.7374,0.7974]';
rn= [0.7913,0.6704,0.6233,0.585]';

bw=	[0.595,0.8835,0.8749,0.3661]';
lw= [0.5837,0.5676,0.4758,0.3288]';
rw= [0.6858,0.6902,0.5093,0.304]';


%%
boxplot([bn,bw,ln,lw,rn,rw],'Labels',{'Bladder(Narrow)','Bladder(Wide)'...
    ,'L.Kidney(Narrow)','L.Kidney(Wide)','R.Kidney(Narrow)','R.Kidney(Wide)'});

%%
[p1,h1] = signrank(bn,bw)
%%
[p,h] = signrank(ln,lw)
%%
[p1,h1] = signrank(rn,rw)

%%
[p1,h1] = ranksum(wn,ww)

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
out = Imap2;
out(Imap2 == 4) = 0;
%imagesc(out(110:290,150:330,200)');
%imagesc(out(260:370,250:360,73)');
imagesc(out(:,:,200)');

axis tight equal off
colormap(map)
caxis([0 4])

%%
imagesc(wM3E1(260:370,250:360,73)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

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
temp = zeros([544 544]);
temp(M1GT(:,:,73) == 1) = 1;
v = [1,1];
contour(temp',v);
axis tight equal off
%%
slice = 130;
imagesc(M4E2(:,:,slice)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
slice = 200;
subplot(1,2,1)
imagesc(GC(110:290,150:340,slice)');
axis tight equal off
caxis([0 4])
colormap(map)

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
Imap = load_raw('C:\Users\yourb\Desktop\new3\ImapM3.raw','*uint8');
%ImapW = load_raw('C:\Users\yourb\Desktop\new3\ImapWM4.raw','*uint8');
GC = load_raw('C:\Users\yourb\Desktop\new3\GC_Mouse3.raw','*uint8');
%GC2 = load_raw('C:\Users\yourb\Desktop\new3\GC2M4.raw','*uint8');
GT = M3GT;
Imap = reshape(Imap,siz); GC = reshape(GC,siz); 
%ImapW = reshape(ImapW,siz); GC2 = reshape(GC2,siz);
%Imap(Imap == 4) = 0;  GT(GT == 4) = 0;  ImapW(ImapW == 4) = 0;
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