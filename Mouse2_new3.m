Xtr = [[pM1E2(pmask1); pM3E2(pmask3)] [pM1E3(pmask1); pM3E3(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E2(pmask2) pM2E3(pmask2) pM2E4(pmask2)];
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
XGTte = [pM2GT(pmask2)];
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
clearvars atlas
atlas  = atlasfunc2(sig1,sig1,K,siz,pmask2,pM1GT,pM3GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood] = ...
    AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask2,siz,30);
JI1= CalcuJI(Imap,pM2GT,K-1);
disp("EM_MAP result")
disp(JI1);
%%
temp = and(Lkidmask,Rkidmask);
temp2 = kidGT - temp;
%%
imagesc(Rkidmask(:,:,200)');
axis tight equal off
%%
blamask = zeros(siz);
temp = Imap == 1; temp2 = zeros(siz);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1; 
temp = bwdist(logical(temp2)) <  power(bwarea(temp(:))/4/pi*3,1/3);
blamask(temp) = 1;
blamask = and(blamask,pmask2); blamask = logical(blamask);

%%
temp = Imap == 2; temp2 = zeros(siz);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1;
Lkidregion =  power(bwarea(temp2(:))/4/pi*3,1/3);


temp = Imap == 3;  temp3 = zeros(siz);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp3(L1.PixelIdxList{idx}) = 1;
Rkidregion =  power(bwarea(temp3(:))/4/pi*3,1/3);

kidregin = (Lkidregion + Rkidregion) /2;
%%
Lkidmask = zeros(siz);
temp = bwdist(logical(temp2)) <  Lkidregion;
Lkidmask(temp) = 1;
Lkidmask = and(Lkidmask,pmask2);
Lkidmask = logical(Lkidmask);

Rkidmask = zeros(siz);
temp = bwdist(logical(temp3)) < Lkidregion;
Rkidmask(temp) = 1;
Rkidmask = and(Rkidmask,pmask2);
Rkidmask = logical(Rkidmask);
%%
temp = or(temp2,temp3);
region =  power(bwarea(temp(:))/4/pi*3,1/3);
temp = bwdist(logical(temp)) < region;
kidmask = and(temp,pmask2);
%%
imagesc(kidmask(:,:,200)');
axis tight equal off
caxis([0 4])
%%
imagesc(Imap(:,:,200)');
axis tight equal off
caxis([0 4])
%%
temp = cutM2GT == 1; temp2 = cutM2GT == 5;
val1 = sum(temp(:)); val2 = sum(temp2(:));
blaratio = val1/ (val1 + val2);

temp = cutM2GT == 2; temp2 = cutM2GT == 6;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Lkidratio = val1/ (val1 + val2);

temp = cutM2GT == 3; temp2 = cutM2GT == 7;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Rkidratio = val1/ (val1 + val2);

%%
Xtebla  = [pM2E2(blamask) pM2E3(blamask) pM2E4(blamask)];
Xtekid = [pM2E2(kidmask) pM2E3(kidmask) pM2E4(kidmask)];

%%
clearvars atlasbla atlaskid
atlasbla   = atlasfunc3(sig1,siz,pmask2,blamask,GMMpro,blaratio,1);
atlaskid = atlasfunc4(sig1,siz,pmask2,kidmask,GMMpro,Lkidratio,Rkidratio,2,3);

%%
GT = zeros(siz);  GT(blamask) = cutM2GT(blamask);
blaGT = zeros(siz); blaGT(blamask) = 3;
blaGT(GT == 1) = 1; blaGT(GT == 5) = 2; 

GT = zeros(siz);  GT(kidmask) = cutM2GT(kidmask);
kidGT = zeros(siz); kidGT(kidmask) = 5;
kidGT(GT == 2) = 1; kidGT(GT == 6) = 2;
kidGT(GT == 3) = 3; kidGT(GT == 7) = 4; 
%%
imagesc(kidGT(:,:,200)');
%%
tmp1 = pM2E2; tmp2 = pM2E3; tmp3 = pM2E4; mask = blaGT;
for k = 1:3
    Sbla.mu(k,1) = mean(tmp1(mask == k));
    Sbla.mu(k,2) = mean(tmp2(mask == k));
    Sbla.mu(k,3) = mean(tmp3(mask == k));
    Sbla.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

mask = kidGT;
for k = 1:5
    Skid.mu(k,1) = mean(tmp1(mask == k));
    Skid.mu(k,2) = mean(tmp2(mask == k));
    Skid.mu(k,3) = mean(tmp3(mask == k));
    Skid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

%%
[Imapbla,~,~,GMMMubla,GMMSigmabla,~,~,~]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,3,blamask,siz,30);
%%
[Imapkid,~,~,GMMMukid,GMMSigmakid,~,~,~]...
    = AtlasGuidedEM_kubo(Xtekid,atlaskid,Skid,5,kidmask,siz,30);

%%
Imap2 = zeros(siz);
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;
Imap2(Imapkid == 1) = 2; Imap2(Imapkid == 2) = 2;
Imap2(Imapkid == 3) = 3; Imap2(Imapkid == 4) = 3;
JI2= CalcuJI(Imap2,pM2GT,K-1);
disp(JI2);

%%
temp = zeros(siz);
temp(Rkidmask) = atlasRkid(:,1);
%%
imagesc(temp(:,:,185)');
axis tight equal
caxis([0 1])
%%
imagesc(Imap2(:,:,180)');
axis tight equal
caxis([0 7])

%%
%In = pM2E2; InGT = blaGT;
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
%mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = pM2E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);

In = pM2E2; InGT = kidGT;
%mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
mu = GMMMukid; sigma = sqrt(GMMSigmakid);

edge =[0 0:0.01:0.9 0.9];
xlim([0 0.9])

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==5),edge,'Normalization','pdf','EdgeAlpha',0.4);

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
plot(edge,y4,'Color',[255 135 0]/255,'LineWidth',2)

mutest = mu(5,1); sigtest = sigma(1,1,5);
y5 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y5,'Color',[255 255 0]/255,'LineWidth',2)
hold off
%%
imagesc(kidGT(:,:,245)');
%%
SE = strel('sphere',1);
temp = zeros(siz);
temp(pM1GT == 1) = 1;
%%
temp2 = imdilate(temp,SE);
%%
imagesc(temp2(:,:,200)');
axis tight equal off