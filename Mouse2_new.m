Xtr = [[pM1E2(pmask1); pM3E2(pmask3)] [pM1E3(pmask1); pM3E3(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E2(pmask2) pM2E3(pmask2) pM2E4(pmask2)];
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
XGTte = [pM2GT(pmask2)];

%%
imagesc(pM2GT(:,:,115)');
axis tight equal off
%%
%initial_value
K=4;
sig1 = 5; %bladder
sig2 = 3; %kidneys
for k = 1:4
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); 
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));

end
clearvars tmp1 tmp2 tmp3

%%
tmp1 = pM2E2; tmp2 = pM2E3; tmp3 = pM2E4; mask = cutGTM2;
for k = 1:7
    S.mu(k,1) = mean(tmp1(mask == k));
    S.mu(k,2) = mean(tmp2(mask == k));
    S.mu(k,3) = mean(tmp3(mask == k));
    S.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3

%%
%Atlas_guided EM2
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask2,pM1GT,pM3GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood] = ...
    AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask2,siz2,30);
JI1= CalcuJI(Imap,pM2GT,K-1);
disp("EM_MAP result")
disp(JI1);

%%
newmaskM2 = zeros(siz2);
for n = 1:3
    temp = Imap == n;
    radi = power(bwarea(temp(:))/4/pi*3,1/3);
    D = bwdist(temp,'euclidean');
    temp = D < radi;
    newmaskM2(temp) = 1;
end
newmaskM2 = and(newmaskM2,pmask2);
newmaskM2 = logical(newmaskM2);
%%
temp = zeros(siz2);
temp(newmaskM2) = pM2GT(newmaskM2);
%%
imagesc(temp(:,:,80)');
axis tight equal off
%%
temp = cutGTM2 == 1; temp2 = cutGTM2 == 5;
val1 = sum(temp(:)); val2 = sum(temp2(:));
blaratio = val1/ (val1 + val2);

temp = cutGTM2 == 2; temp2 = cutGTM2 == 6;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Lkidration = val1/ (val1 + val2);

temp = cutGTM2 == 3; temp2 = cutGTM2 == 7;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Rkidration = val1/ (val1 + val2);
%%
imagesc(Imap2(:,:,200)');
caxis([0 4])
colormap(map)
axis tight equal off
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
%%
temp = zeros(siz2);
temp(newmaskM2) = atlasnew(:,2);
%%
imagesc(cutGTM2(:,:,200)');
caxis([0 7])
axis tight equal off
%%
%Atlas_guided EM
clearvars atlasnew
atlasnew  = atlasfunc3(sig1,sig2,siz2,pmask2,newmaskM2,GMMpro,blaratio,Lkidration,Rkidration);
Snew.mu = GMMMu; Snew.Sigma = GMMSigma;
Xtenew = [pM2E2(newmaskM2) pM2E3(newmaskM2) pM2E4(newmaskM2)];
%%
[Imap2,~,PP2,GMMMu2,GMMSigma2,GMMpro2,~,likelihood2] ...
    = AtlasGuidedEM_kubo(Xtenew,atlasnew,S,7,newmaskM2,siz2,30);
JI2= CalcuJI(Imap2,pM2GT,K-1);
disp("EM_MAP result")
disp(JI2);
%%
Imap3 = Imap2;
Imap3(Imap2 == 5) = 1;
Imap3(Imap2 == 6) = 2;
Imap3(Imap2 == 7) = 3;
JI2= CalcuJI(Imap3,pM2GT,K-1);
disp("EM_MAP result")
disp(JI2);
%%
imagesc(temp(:,:,80)');
axis tight equal off
caxis([0 7])
%%
imagesc(pM2E2(:,:,200)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
In = pM2E2; InGT = cutGTM2;
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
In = pM2E2; InGT = cutGTM2;
edge =[0 0:0.01:1.2 1.2];
hold on
histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);
histogram(In(temp ==4),edge,'EdgeAlpha',0.4);

histogram(In(InGT ==5),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==6),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==7),edge,'EdgeAlpha',0.4);
%%
In = pM2E2; InGT = pM2GT;
%mu = SS.mu; sigma = sqrt(SS.Sigma);
mu = GMMMu; sigma = sqrt(GMMSigma);

edge =[0 0:0.01:0.9 0.9];
xlim([0 0.9])

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4,'FaceColor',[51 102 255]/255);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4,'FaceColor',[255 135 0]/255);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4,'FaceColor',[255 255 0]/255);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4,'FaceColor',[142 0 204]/255);

mutest = mu(1,1); sigtest = sigma(1,1,1);
y1 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest = mu(2,1); sigtest =  sigma(1,1,2);
y2 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest = mu(3,1); sigtest = sigma(1,1,3);
y3 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

mutest = mu(4,1); sigtest =  sigma(1,1,4);
y4 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y4,'Color',[142 0 204]/255,'LineWidth',2)