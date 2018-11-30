%train_mouse2_mouse3 test_mouse1
Xtr = [[pM2E2(pmask2); pM3E2(pmask3)] [pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E2(pmask1) pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];
XGTte = [pM1GT(pmask1)];
%%
%initial_value
K = 4;
for k = 1:K
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
 clearvars tmp1 tmp2 tmp3
 
%%
tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4; mask = cutGTM1;
for k = 1:3
    S.mu(k,1) = mean(tmp1(mask == k));
    S.mu(k,2) = mean(tmp2(mask == k));
    S.mu(k,3) = mean(tmp3(mask == k));
    S.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
    S.mu(k,:) = S.mu(k,:) - 5*sqrt(diag(S.Sigma(:,:,k)))';
    S.Sigma(:,:,k) = (sqrt(S.Sigma(:,:,k)).*2).^2;
end
clearvars tmp1 tmp2 tmp3
%%
%Atlas_guided EM
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask1,pM2GT,pM3GT);
%%
[Imap,~,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask1,siz2,30);
JI1= CalcuJI(Imap,pM1GT,K-1);
disp("EM_MAP result")
disp(JI1);
%%
temp = Imap ==1;
BWM1 = zeros(siz2);
L1 = bwconncomp(temp);
numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels);
BWM1(L1.PixelIdxList{idx}) = 1;
%%
imagesc(BWM1(:,:,80)');
%%
blamask = zeros(siz2);

temp = Imap == 1; temp = logical(temp); temp2 = zeros(siz2);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1;
D = bwdist(logical(temp2),'euclidean');
temp = D <  power(bwarea(temp(:))/4/pi*3,1/3);
blamask(temp) = 1;
blamask = and(blamask,pmask1); blamask = logical(blamask);
%%
Rkidmask = zeros(siz2);

temp = Imap == 2; temp = logical(temp); temp2 = zeros(siz2);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1;
D = bwdist(logical(temp2),'euclidean');
temp = D <  power(bwarea(temp(:))/4/pi*3,1/3);
Rkidmask(temp) = 1;
Rkidmask = and(Rkidmask,pmask1); Rkidmask = logical(Rkidmask);

%%
Lkidmask = zeros(siz2);

temp = Imap == 3; temp = logical(temp); temp2 = zeros(siz2);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1;
D = bwdist(logical(temp2),'euclidean');
temp = D <  power(bwarea(temp(:))/4/pi*3,1/3);
Lkidmask(temp) = 1;
Lkidmask = and(Lkidmask,pmask1); Lkidmask = logical(Lkidmask);
%%
imagesc(GT1(:,:,105)');
%%
GT1 = zeros(siz2);
GT1(blamask) = pM1GT(blamask);
%%
temp = cutGTM1 == 1; temp2 = cutGTM1 == 5;
val1 = sum(temp(:)); val2 = sum(temp2(:));
blaratio = val1/ (val1 + val2);

temp = cutGTM1 == 2; temp2 = cutGTM1 == 6;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Lkidration = val1/ (val1 + val2);

temp = cutGTM1 == 3; temp2 = cutGTM1 == 7;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Rkidration = val1/ (val1 + val2);
%%
%Atlas_guided EM
clearvars atlasnew
atlasnew  = atlasfunc3(sig1,sig2,siz2,pmask1,blamask,GMMpro,blaratio,Lkidration,Rkidration);
Snew.mu = GMMMu; Snew.Sigma = GMMSigma;
Xtenew = [pM1E2(blamask) pM1E3(blamask) pM1E4(blamask)];
%%
[Imap2,~,PP2,GMMMu2,GMMSigma2,GMMpro2,~,likelihood2] ...
    = AtlasGuidedEM_kubo(Xtenew,atlasnew,S,7,blamask,siz2,30);
JI2= CalcuJI(Imap2,pM1GT,K-1);
disp("EM_MAP result")

disp(JI2);
%%
Imap3 = Imap2;
Imap3(Imap2 == 5) = 1;
Imap3(Imap2 == 6) = 2;
Imap3(Imap2 == 7) = 3;
JI2= CalcuJI(Imap3,pM1GT,K-1);
disp("EM_MAP result")
disp(JI2);
%%
temp = zeros(siz2);
temp(blamask) = atlasnew(:,2);
%%
imagesc(blamask(:,:,75)');
axis tight equal off
%%
tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4; mask = cutGTM1;
for k = 1:7
    S.mu(k,1) = mean(tmp1(mask == k));
    S.mu(k,2) = mean(tmp2(mask == k));
    S.mu(k,3) = mean(tmp3(mask == k));
    
    S.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
imagesc(pM1E2(:,:,70)');
caxis([0 0.7])
colormap(gray)
axis tight equal off
%%
In = pM1E2; InGT = cutGTM1;
mu = S.mu; sigma = sqrt(S.Sigma);
%mu = GMMMu; sigma = sqrt(GMMSigma);
%mu = GMMMu2; sigma = sqrt(GMMSigma2);
edge =[0 0:0.01:2.79 2.8];
xlim([0 2.8])

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
In = pM1E2; InGT = cutGTM1;
mu = S.mu; sigma = sqrt(S.Sigma);
%mu = GMMMu; sigma = sqrt(GMMSigma);
%mu = GMMMu2; sigma = sqrt(GMMSigma2);
edge =[0 0:0.01:2.8 2.8];
xlim([0 2.8])

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
In = pM1E2; InGT = GT1;
edge =[0 0:0.01:2.7 2.7];
hold on
histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'EdgeAlpha',0.4);
%%
histogram(In(InGT ==5),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==6),edge,'EdgeAlpha',0.4);