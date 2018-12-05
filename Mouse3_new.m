%train_mouse1 mouse2 mouse4 test_mouse3
Xtr = [[M1E2(mask1); M2E2(mask2) M4E2(mask4)] [M1E3(mask1); M2E3(mask2) M4E3(mask4)]...
      [M1E4(mask1); M2E4(mask2)] M4E4(mask4)];
Xte = [M3E2(mask3) M3E3(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2) M4GT(mask4)];
XGTte = M3GT(mask3);
%%
%initial_value
K=4;
sig1 = 5; %bladder
sig2 = 3; %kidneys

for k = 1:K
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); 
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
%initial_value test
for k = 1:K
    tmp1 = Xte(:,1); tmp2 = Xte(:,2); tmp3 = Xte(:,3);
    STe.mu(k,1) = mean(tmp1(XGTte == k));
    STe.mu(k,2) = mean(tmp2(XGTte == k));
    STe.mu(k,3) = mean(tmp3(XGTte == k));
    STe.Sigma(:,:,k) = cov(([tmp1(XGTte == k),tmp2(XGTte == k),tmp3(XGTte == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
%Atlas_guided EM
atlas  = atlasfunc2(sig1,sig2,K,siz2,mask3,M1GT,M2GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,lilelihood] ...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,mask3,siz2,30);
JI1= CalcuJI(Imap,M3GT,K-1);
disp("EM_MAP result")
disp(JI1);
%%
imagesc(Imap(:,:,155)');
axis tight equal off
caxis([0 4])
%%
newmaskM3 = zeros(siz2);
for n = 1:3
    temp = Imap == n;
    radi = power(bwarea(temp(:))/4/pi*3,1/3);
    D = bwdist(temp,'euclidean');
    temp = D < radi;
    newmaskM3(temp) = 1;
end
newmaskM3 = and(newmaskM3,mask3);
newmaskM3 = logical(newmaskM3);
%%
GT3 = zeros(siz2);
GT3(newmaskM3) = M3GT(newmaskM3);
%%
tmp1 = pM3E2; tmp2 = pM3E3; tmp3 = pM3E4; mask = cutGTM3;
for k = 1:7
    S.mu(k,1) = mean(tmp1(mask == k));
    S.mu(k,2) = mean(tmp2(mask == k));
    S.mu(k,3) = mean(tmp3(mask == k));
    S.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
temp = cutGTM3 == 1; temp2 = cutGTM3 == 5;
val1 = sum(temp(:)); val2 = sum(temp2(:));
blaratio = val1/ (val1 + val2);

temp = cutGTM3 == 2; temp2 = cutGTM3 == 6;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Lkidration = val1/ (val1 + val2);

temp = cutGTM3 == 3; temp2 = cutGTM3 == 7;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Rkidration = val1/ (val1 + val2);

%%
%Atlas_guided EM
clearvars atlasnew
atlasnew  = atlasfunc3(sig1,sig2,siz2,mask3,newmaskM3,GMMpro,blaratio,Lkidration,Rkidration);
Snew.mu = GMMMu; Snew.Sigma = GMMSigma;
Xtenew = [pM3E2(newmaskM3) pM3E3(newmaskM3) pM3E4(newmaskM3)];
%%
[Imap2,~,PP2,GMMMu2,GMMSigma2,GMMpro2,~,likelihood2] ...
    = AtlasGuidedEM_kubo(Xtenew,atlasnew,S,7,newmaskM3,siz2,30);
JI2= CalcuJI(Imap2,M3GT,K-1);
disp("EM_MAP result")
disp(JI2);
%%
Imap3 = Imap2;
Imap3(Imap2 == 5) = 1;
Imap3(Imap2 == 6) = 2;
Imap3(Imap2 == 7) = 3;
JI2= CalcuJI(Imap3,M3GT,K-1);
disp("EM_MAP result")
disp(JI2);
%%
imagesc(pM3E2(:,:,200)');
axis tight equal off
caxis([0 0.7]);
colormap(gray)
%%
In = pM3E2; InGT = cutGTM3;
%mu = S.mu; sigma = sqrt(S.Sigma);
%mu = GMMMu; sigma = sqrt(GMMSigma);
mu = GMMMu2; sigma = sqrt(GMMSigma2);
edge =[0 0:0.01:1.2 1.2];
xlim([0 1.2])

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);

histogram(In(InGT ==5),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==6),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==7),edge,'Normalization','pdf','EdgeAlpha',0.4);

mutest = mu(1,1); sigtest = sigma(1,1,1);
y1 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest = mu(2,1); sigtest =  sigma(1,1,2);
y2 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest = mu(3,1); sigtest = sigma(1,1,3);
y3 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

mutest = mu(4,1); sigtest = sigma(1,1,4);
y4 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y4,'Color',[142 0 204]/255,'LineWidth',2)

mutest = mu(5,1); sigtest = sigma(1,1,5);
y5 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y5,'Color',[35 167 22]/255,'LineWidth',2)

mutest = mu(6,1); sigtest = sigma(1,1,6);
y6 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y6,'Color',[15 82 188]/255,'LineWidth',2)

mutest = mu(7,1); sigtest = sigma(1,1,7);
y7 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y7,'Color',[204 0 0]/255,'LineWidth',2)

%%
In = pM3E2; InGT = cutGTM3;
edge =[0 0:0.01:2.5 2.5];
hold on
histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==5),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==6),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==7),edge,'EdgeAlpha',0.4);
%%
In = pM3E2; InGT = GT3;
edge =[0 0:0.01:2.0 2.0];
hold on
histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'EdgeAlpha',0.4);
