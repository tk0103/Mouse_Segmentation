%train_mouse2_mouse3 test_mouse1
Xtr = [[pM2E2(pmask2); pM3E2(pmask3)] [pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E2(pmask1) pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];
XGTte = [pM1GT(pmask1)];
mask = pmask1;

%%
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
%Atlas_guided EM
atlas  = atlasfunc_old(sig1,sig2,K,siz,mask,pM2GT,pM3GT);
%%
[Imap,~,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,mask,siz,30);
%%
JI1= CalcuJI(Imap,pM1GT,K-1);
disp("EM_MAP result")
disp(JI1);
%%
imagesc(pM1GT(:,:,70)');
axis tight equal
%%
blamask = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(tmp) = 1;
blamask = logical(and(blamask,pmask1)); 

Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L1 = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(tmp) = 1;
Lkidmask = logical(and(Lkidmask,pmask1));

Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L1 = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask(tmp) = 1;
Rkidmask = logical(and(Rkidmask,pmask1));
LRAND = and(Rkidmask,Lkidmask);
%%
[XX,YY,ZZ] = meshgrid(1:siz(1),1:siz(2),1:siz(3));
tmp = bwdist(logical(not(Lmaxcomp)));
[~,I] = max(tmp(:)); [y,x,z] = ind2sub(siz,I);
RLkid = sqrt((XX - x).^2 + (YY -y).^2 + (ZZ - z).^2);

tmp = bwdist(logical(not(Rmaxcomp)));
[~,I] = max(tmp(:)); [y,x,z] = ind2sub(siz,I);
RRkid = sqrt((XX - x).^2 + (YY -y).^2 + (ZZ - z).^2);

tmp = zeros(siz);
tmp(RLkid > RRkid) = 1;
tmp2 = and(tmp,LRAND);
Rkidmask = logical(Rkidmask - LRAND + tmp2); 

tmp(RRkid > RLkid) = 1;
tmp2 = and(tmp,LRAND);
Lkidmask = logical(Lkidmask - LRAND + tmp2);
%%
temp = cutM1GT == 1; temp2 = cutM1GT == 5;
val1 = sum(temp(:)); val2 = sum(temp2(:));
blaratio = val1/ (val1 + val2);

temp = cutM1GT == 2; temp2 = cutM1GT == 6;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Lkidration = val1/ (val1 + val2);

temp = cutM1GT == 3; temp2 = cutM1GT == 7;
val1 = sum(temp(:)); val2 = sum(temp2(:));
Rkidration = val1/ (val1 + val2);

%%
XGTtr = [blaGT(blamask); pM3GT(blamask)];
%%
K = 3;
for k = 1:K
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
 clearvars tmp1 tmp2 tmp3
%%
Xtebla = [pM1E2(blamask) pM1E3(blamask) pM1E4(blamask)];
XteLkid = [pM1E2(Lkidmask) pM1E3(Lkidmask) pM1E4(Lkidmask)];
XteRkid = [pM1E2(Rkidmask) pM1E3(Rkidmask) pM1E4(Rkidmask)];
%%
atlasbla   = atlasfunc2(sig2,siz,pmask1,blamask,GMMpro,0.8,1);
atlasLkid  = atlasfunc2(sig2,siz,pmask1,Lkidmask,GMMpro,0.8,2);
atlasRkid  = atlasfunc2(sig2,siz,pmask1,Rkidmask,GMMpro,0.8,3);

%%
GT = zeros(siz);  GT(blamask) = cutM1GT(blamask);
blaGT = zeros(siz); blaGT(blamask) = 3;
blaGT(GT == 1) = 1; blaGT(GT == 5) = 2; 

GT = zeros(siz);  GT(Lkidmask) = cutM1GT(Lkidmask);
LkidGT = zeros(siz); LkidGT(Lkidmask) = 3;
LkidGT(GT == 2) = 1; LkidGT(GT == 6) = 2;

GT = zeros(siz);  GT(Rkidmask) = cutM1GT(Rkidmask);
RkidGT = zeros(siz); RkidGT(Rkidmask) = 3;
RkidGT(GT == 3) = 1; RkidGT(GT == 7) = 2; 
%%
imagesc(RkidGT(:,:,195)');
axis tight equal off
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

SRkid.mu(1,:) = GMMMu(3,:) -0.25*sqrt(diag(GMMSigma(:,:,3)))';
SRkid.mu(2,:) = GMMMu(3,:) +0.25*sqrt(diag(GMMSigma(:,:,3)))';
SRkid.mu(3,:) = GMMMu(4,:); 

SRkid.Sigma(:,:,1) =  (sqrt(GMMSigma(:,:,3))./4).^2;
SRkid.Sigma(:,:,2) =  (sqrt(GMMSigma(:,:,3))./4).^2;
SRkid.Sigma(:,:,3) =  (sqrt(GMMSigma(:,:,4))).^2;
%%
tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4; mask = blaGT;
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
[Imapbla,~,~,GMMMubla,GMMSigmabla,GMMprobla,Featbla,likelihoodbla]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,3,blamask,siz,30);
%%
p_l = atlasbla;
r = pM1E2 <  GMMMubla(2,1); 
r = r(blamask);
p_l(r,2) = 0;

PP = Featbla.*p_l;
PP = bsxfun(@rdivide,PP,sum(PP,2));

[~,L] = max(PP,[],2);
Imapbla = zeros(siz);
Imapbla(blamask) = L;
%%
[ImapLkid,~,~,GMMMuLkid,GMMSigmaLkid,~,~,likelihoodLkid]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,3,Lkidmask,siz,30);
%%
[ImapRkid,~,~,GMMMuRkid,GMMSigmaRkid,~,~,likelihoodRkid]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,3,Rkidmask,siz,30);
%%
Imap2 = zeros(siz);
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;

JI2= CalcuJI(Imap2,pM1GT,K-1);
disp(JI2);
%%
imagesc(Imap2(:,:,200)');
axis tight equal off
caxis([0 4])
%colormap(gray)
%%
imagesc(pM2E2(:,:,220)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
Imap2 = zeros(siz);
Imap2(Imapbla == 1) = 1;
Imap2(Imapbla == 2) = 1;

Imap2(ImapLkid == 1) = 2;
Imap2(ImapLkid == 2) = 2;

Imap2(ImapRkid == 1) = 3;
Imap2(ImapRkid == 2) = 3;

JI2= CalcuJI(Imap2,pM1GT,K-1);
disp(JI2);
%%
In = pM1E2; InGT = blaGT;
mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
%mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = pM1E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);

%In = pM1E2; InGT = RkidGT;
%mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
%mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.01:1.7 1.7];
xlim([0 1.7])

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
edge =[0 0:0.01:3.0 3.0];
 In = pM1E2; InGT = pM1GT;
hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);

%%
In = pM1E2; InGT = pM1GT;
%mu = SS.mu; sigma = sqrt(SS.Sigma);
mu = GMMMu; sigma = sqrt(GMMSigma);

edge =[0 0:0.01:2.7 2.7];
xlim([0 2.7])

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