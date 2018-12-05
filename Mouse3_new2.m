Xtr = [[pM1E2(pmask1); pM2E2(pmask2)] [pM1E3(pmask1); pM2E3(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)] ];
Xte = [pM3E2(pmask3) pM3E3(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];
XGTte = pM3GT(pmask3);
mask = pmask3; mask = logical(mask);
%%
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
atlas  = atlasfunc2(sig1,sig2,K,siz,mask,pM1GT,pM2GT);

%%
[Imap,~,PP,GMMMu,GMMSigma,GMMpro,Feat,lilelihood] ...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,mask,siz,30);
JI1= CalcuJI(Imap,pM3GT,K-1);
disp("EM_MAP result")
disp(JI1);
%%
imagesc(Imap(:,:,70)');
axis tight equal off
caxis([0 4])
%%
blamask = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
temp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(temp) = 1;
blamask = logical(and(blamask,mask)); 

Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L1 = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
temp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(temp) = 1;
Lkidmask = logical(and(Lkidmask,mask));

Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L1 = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
temp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask(temp) = 1;
Rkidmask = logical(and(Rkidmask,mask));

temp = zeros(siz);
temp(bwdist(logical(Rmaxcomp)) > bwdist(logical(Rmaxcomp))) = 1;
temp2 = and(temp,and(Lkidmask,Rkidmask));
Rkidmask = logical(Rkidmask - temp2); 

temp(bwdist(logical(Rmaxcomp)) < bwdist(logical(Lmaxcomp))) = 1;
temp2 = and(temp,and(Lkidmask,Rkidmask));
Lkidmask = logical(Lkidmask - temp2);
%%
temp = cutM3GT == 1; temp2 = cutM3GT == 5;
val1 = sum(temp(:)); val2 = sum(temp2(:));
blaratio = val1/ (val1 + val2);

temp = cutM3GT == 2; temp2 = cutM3GT == 6;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Lkidration = val1/ (val1 + val2);

temp = cutM3GT == 3; temp2 = cutM3GT == 7;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Rkidration = val1/ (val1 + val2);

%%
Xtebla  = [pM3E2(blamask)  pM3E3(blamask) pM3E4(blamask)];
XteLkid = [pM3E2(Lkidmask) pM3E3(Lkidmask) pM3E4(Lkidmask)];
XteRkid = [pM3E2(Rkidmask) pM3E3(Rkidmask) pM3E4(Rkidmask)];
%%
clearvars atlasbla atlasLkid atlasRkid
atlasbla   = atlasfunc3(sig1,siz,mask,blamask,GMMpro,0.8,1);
atlasLkid  = atlasfunc3(sig1,siz,mask,Lkidmask,GMMpro,0.8,2);
atlasRkid  = atlasfunc3(sig1,siz,mask,Rkidmask,GMMpro,0.8,3);

%%
atlastemp1 = zeros(siz); 
atlastemp1(mask) = GMMpro(:,1);
atlastemp2 = zeros(siz); 
atlastemp2(mask) = GMMpro(:,2);
atlastemp3 = zeros(siz); 
atlastemp3(mask) = GMMpro(:,3);

atlastemp4 = zeros(siz); 
atlastemp4(mask) = GMMpro(:,4);
%%
atlastemp1 = squeeze(atlastemp1(:,:,200));
atlastemp2 = squeeze(atlastemp2(:,:,200));
atlastemp3 = squeeze(atlastemp3(:,:,200));
atlastemp4 = squeeze(atlastemp4(:,:,200));
newmask = mask(:,:,200);
%%
atat = zeros(49904,3);
%atat(:,1) = atlastemp1(newmask);
atat(:,1) = atlastemp2(newmask);
atat(:,2) = atlastemp3(newmask);
atat(:,3) = atlastemp4(newmask);
%%
Xte = [pM3E2(newmask) pM3E3(newmask) pM3E4(newmask)];
%%
temp = zeros(siz);
temp(pmask3) = atlasRkid(:,3);

%%
r = atlasRkid < 0;
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,lilelihood] ...
    = AtlasGuidedEM_kubo(Xte,atat,SS,3,newmask,siz,30);
%%
imagesc(temp(:,:,200)');
%%
atlasRkid  = atlasfunc3(sig1,siz,mask,Rkidmask,GMMpro,0.8,3);
%%
GT = zeros(siz);  GT(blamask) = cutM3GT(blamask);
blaGT = zeros(siz); blaGT(blamask) = 3;
blaGT(GT == 1) = 1; blaGT(GT == 5) = 2; 

GT = zeros(siz);  GT(Lkidmask) = cutM3GT(Lkidmask);
LkidGT = zeros(siz); LkidGT(Lkidmask) = 3;
LkidGT(GT == 2) = 1; LkidGT(GT == 6) = 2;

GT = zeros(siz);  GT(Rkidmask) = cutM3GT(Rkidmask);
RkidGT = zeros(siz); RkidGT(Rkidmask) = 3;
RkidGT(GT == 3) = 1; RkidGT(GT == 7) = 2; 
%%
imagesc(blaGT(:,:,55)');
%%
Sbla.mu(1,:) = GMMMu(1,:) +0.5*sqrt(diag(GMMSigma(:,:,1)))';
Sbla.mu(2,:) = GMMMu(1,:) -0.5*sqrt(diag(GMMSigma(:,:,1)))';
Sbla.mu(3,:) = GMMMu(4,:); 

Sbla.Sigma(:,:,1) =  (sqrt(GMMSigma(:,:,1))./4).^2;
Sbla.Sigma(:,:,2) =  (sqrt(GMMSigma(:,:,1))./4).^2;
Sbla.Sigma(:,:,3) =  GMMSigma(:,:,4);

SLkid.mu(1,:) = GMMMu(2,:) -0.25*sqrt(diag(GMMSigma(:,:,2)))';
SLkid.mu(2,:) = GMMMu(2,:) +0.25*sqrt(diag(GMMSigma(:,:,2)))';
SLkid.mu(3,:) = GMMMu(4,:); 

SLkid.Sigma(:,:,1) =  (sqrt(GMMSigma(:,:,2))./4).^2;
SLkid.Sigma(:,:,2) =  (sqrt(GMMSigma(:,:,2))./4).^2;
SLkid.Sigma(:,:,3) =  (sqrt(GMMSigma(:,:,4))).^2;
%%
SRkid.mu(1,:) = GMMMu(3,:) -0.25*sqrt(diag(GMMSigma(:,:,3)))';
SRkid.mu(2,:) = GMMMu(3,:) +0.25*sqrt(diag(GMMSigma(:,:,3)))';
SRkid.mu(3,:) = GMMMu(4,:); 

SRkid.Sigma(:,:,1) =  (sqrt(GMMSigma(:,:,3))./4).^2;
SRkid.Sigma(:,:,2) =  (sqrt(GMMSigma(:,:,3))./4).^2;
SRkid.Sigma(:,:,3) =  (sqrt(GMMSigma(:,:,4))).^2;
%%
tmp1 = pM3E2; tmp2 = pM3E3; tmp3 = pM3E4; mask = blaGT;
for k = 1:3
    Sbla.mu(k,1) = mean(tmp1(mask == k));
    Sbla.mu(k,2) = mean(tmp2(mask == k));
    Sbla.mu(k,3) = mean(tmp3(mask == k));
    Sbla.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

mask = LkidGT;
for k = 1:3
    SLkid.mu(k,1) = mean(tmp1(mask == k));
    SLkid.mu(k,2) = mean(tmp2(mask == k));
    SLkid.mu(k,3) = mean(tmp3(mask == k));
    SLkid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

mask = RkidGT;
for k = 1:3
    SRkid.mu(k,1) = mean(tmp1(mask == k));
    SRkid.mu(k,2) = mean(tmp2(mask == k));
    SRkid.mu(k,3) = mean(tmp3(mask == k));
    SRkid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3

%%
[Imapbla,~,~,GMMMubla,GMMSigmabla,~,Featbla,likelihoodbla]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,3,blamask,siz,30);
%%
[ImapLkid,~,~,GMMMuLkid,GMMSigmaLkid,~,~,likelihoodLkid]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,3,Lkidmask,siz,30);
%%
[ImapRkid,~,~,GMMMuRkid,GMMSigmaRkid,~,~,likelihoodRkid]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,3,Rkidmask,siz,30);
%%
temp = zeros(siz);
temp(Rkidmask) = atlasRkid(:,3);
%%
imagesc(temp(:,:,230)');
axis tight equal
caxis([0 1])
%%
p_l = atlasbla;
r = pM3E2 <  GMMMubla(2,1);  r = r(blamask);
p_l(r,2) = 0;
PP = (Featbla.*p_l);
PP = bsxfun(@rdivide,PP,sum(PP,2));

[~,L] = max(PP,[],2);
Imapbla = zeros(siz);
Imapbla(blamask) = L;
%%
Imap2 = zeros(siz);
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
JI2= CalcuJI(Imap2,pM3GT,K-1);

disp(JI2);

%%
imagesc(Imap2(:,:,220)');
axis tight equal off
caxis([0 4])
%%
imagesc(pM3E2(:,:,80)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
%In = pM3E2; InGT = blaGT;
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
%mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = pM3E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);


In = pM3E2; InGT = RkidGT;
mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.01:0.9 0.9];
xlim([0 0.9])

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);

mutest = mu(1,1); sigtest = sigma(1,1,1);
y1 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

mutest = mu(2,1); sigtest =  sigma(1,1,2);
y2 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

mutest = mu(3,1); sigtest = sigma(1,1,3);
y3 = pdf('Normal',edge,mutest,sigtest);
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)
hold off
%%
edge =[0 0:0.01:1.5 1.5];
 In = pM3E2; InGT = RkidGT;
hold on
histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);

%%
SE = strel('sphere',2);
temp = zeros(siz);
temp(pM1GT ==1) = 1;
%%
temp = pM1GT;
temp(temp2 == 1) = 1;
%%
temp2 = imdilate(temp,SE);
%%
pM1GT = temp;
%%
In = pM3E2; InGT = pM3GT;
mu = SS.mu; sigma = sqrt(SS.Sigma);
mu = GMMMu; sigma = sqrt(GMMSigma);

edge =[0 0:0.01:0.7 0.7];
xlim([0 0.7])

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