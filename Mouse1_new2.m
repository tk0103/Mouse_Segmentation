%train_mouse2_mouse3 test_mouse1
Xtr = [[pM2E2(pmask2); pM3E2(pmask3)] [pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E2(pmask1) pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];
XGTte = [pM1GT(pmask1)];
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
atlas  = atlasfunc2(sig1,sig2,K,siz,pmask1,pM2GT,pM3GT);
%%
[Imap,~,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask1,siz,30);
%%
JI1= CalcuJI(Imap,temp,K-1);
disp("EM_MAP result")
disp(JI1);

%%
blamask = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
temp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(temp) = 1;
blamask = logical(and(blamask,pmask1)); 
%%
Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L1 = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
temp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(temp) = 1;
Lkidmask = logical(and(Lkidmask,pmask1));

Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L1 = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
temp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask(temp) = 1;
Rkidmask = logical(and(Rkidmask,pmask1));
%%
temp = and(Lkidmask,Rkidmask);
%%
temp2 = zeros(siz);
temp2(bwdist(logical(Rmaxcomp)) > bwdist(logical(Lmaxcomp))) = 1;
%%

imagesc(temp2(:,:,200)');


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
atlasbla   = atlasfunc3(sig2,siz,pmask1,blamask,GMMpro,blaratio,1);
atlasLkid  = atlasfunc3(sig2,siz,pmask1,Lkidmask,GMMpro,Lkidration,2);
atlasRkid  = atlasfunc3(sig2,siz,pmask1,Rkidmask,GMMpro,Rkidration,3);

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
[Imapbla,~,~,GMMMubla,GMMSigmabla,~,~,~]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,3,blamask,siz,30);

[ImapLkid,~,~,GMMMuLkid,GMMSigmaLkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,3,Lkidmask,siz,30);

[ImapRkid,~,~,GMMMuRkid,GMMSigmaRkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,3,Rkidmask,siz,30);
%%

Imap2 = zeros(siz);
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;

JI2= CalcuJI(Imap2,temp,K-1);
disp(JI2);
%%
imagesc(Imap2(:,:,70)');
axis tight equal off
caxis([0 4])
%colormap(gray)
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
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = pM1E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);

%In = pM1E2; InGT = RkidGT;
% mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);

%mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.01:2.5 2.5];
xlim([0 2.5])

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

%%
SE = strel('sphere',1);
temp = zeros(siz);
temp(pM1GT ==1) = 1;
%%
temp = pM1GT;
temp(Lmaxcomp == 1) = 1;
%%
Lmaxcomp = imdilate(temp,SE);