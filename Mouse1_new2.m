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
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask1,pM2GT,pM3GT);
%%
[Imap,~,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask1,siz2,30);
JI1= CalcuJI(Imap,pM1GT,K-1);
disp("EM_MAP result")
disp(JI1);

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
Lkidmask = zeros(siz2);

temp = Imap == 2; temp = logical(temp); temp2 = zeros(siz2);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1;
D = bwdist(logical(temp2),'euclidean');
temp = D <  power(bwarea(temp(:))/4/pi*3,1/3);
Lkidmask(temp) = 1;
Lkidmask = and(Lkidmask,pmask1); Lkidmask = logical(Lkidmask);
%%
Rkidmask = zeros(siz2);

temp = Imap == 3; temp = logical(temp); temp2 = zeros(siz2);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1;
D = bwdist(logical(temp2),'euclidean');
temp = D <  power(bwarea(temp(:))/4/pi*3,1/3);
Rkidmask(temp) = 1;
Rkidmask = and(Rkidmask,pmask1); Rkidmask = logical(Rkidmask);


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
Xtebla = [pM1E2(blamask) pM1E3(blamask) pM1E4(blamask)];
XteLkid = [pM1E2(Lkidmask) pM1E3(Lkidmask) pM1E4(Lkidmask)];
XteRkid = [pM1E2(Rkidmask) pM1E3(Rkidmask) pM1E4(Rkidmask)];
%%
atlasbla   = atlasfunc3(sig2,siz2,pmask1,blamask,GMMpro,blaratio,1);
atlasLkid  = atlasfunc3(sig2,siz2,pmask1,Lkidmask,GMMpro,Lkidration,2);
atlasRkid  = atlasfunc3(sig2,siz2,pmask1,Rkidmask,GMMpro,Rkidration,3);
%%
Sbla.mu(1,:) = GMMMu(1,:); Sbla.mu(2,:) = GMMMu(1,:); Sbla.mu(3,:) = GMMMu(4,:);
Sbla.Sigma(:,:,1) = GMMSigma(:,:,1);
Sbla.Sigma(:,:,2) = GMMSigma(:,:,1);
Sbla.Sigma(:,:,3) = GMMSigma(:,:,4);

SLkid.mu(1,:) = GMMMu(2,:); SLkid.mu(2,:) = GMMMu(2,:); SLkid.mu(3,:) = GMMMu(4,:);
SLkid.Sigma(:,:,1) = GMMSigma(:,:,2);
SLkid.Sigma(:,:,2) = GMMSigma(:,:,2);
SLkid.Sigma(:,:,3) = GMMSigma(:,:,4);

SRkid.mu(1,:) = GMMMu(3,:); SRkid.mu(2,:) = GMMMu(3,:); SRkid.mu(3,:) = GMMMu(4,:);
SRkid.Sigma(:,:,1) = GMMSigma(:,:,3);
SRkid.Sigma(:,:,2) = GMMSigma(:,:,3);
SRkid.Sigma(:,:,3) = GMMSigma(:,:,4);
%%
blaGT = zeros(siz2);
blaGT(blamask) = cutGTM1(blamask);
blaGT(blaGT == 5) = 2; blaGT(blaGT == 4) = 3;

LkidGT = zeros(siz2);
LkidGT(Lkidmask) = cutGTM1(Lkidmask);
LkidGT(LkidGT == 2) = 1; LkidGT(LkidGT == 6) = 2; LkidGT(LkidGT == 4) = 3;

RkidGT = zeros(siz2);
RkidGT(Rkidmask) = cutGTM1(Rkidmask);
RkidGT(RkidGT == 3) = 1; RkidGT(RkidGT == 7) = 2; RkidGT(RkidGT == 4) = 3;

%%
tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4; mask = blaGT;
for k = 1:3
    Sbla.mu(k,1) = mean(tmp1(mask == k));
    Sbla.mu(k,2) = mean(tmp2(mask == k));
    Sbla.mu(k,3) = mean(tmp3(mask == k));
    Sbla.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4; mask = LkidGT;
for k = 1:3
    SLkid.mu(k,1) = mean(tmp1(mask == k));
    SLkid.mu(k,2) = mean(tmp2(mask == k));
    SLkid.mu(k,3) = mean(tmp3(mask == k));
    SLkid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4; mask = RkidGT;
for k = 1:3
    SRkid.mu(k,1) = mean(tmp1(mask == k));
    SRkid.mu(k,2) = mean(tmp2(mask == k));
    SRkid.mu(k,3) = mean(tmp3(mask == k));
    SRkid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
[Imapbla,~,~,GMMMubla,GMMSigmabla,~,~,~]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,3,blamask,siz2,30);

[ImapLkid,~,~,GMMMuLkid,GMMSigmaLkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,3,Lkidmask,siz2,30);

[ImapRkid,~,~,GMMMuRkid,GMMSigmaRkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,3,Rkidmask,siz2,30);
%%
temp = zeros(siz2);
temp(Rkidmask) = atlasRkid(:,3);
%%
imagesc(ImapRkid(:,:,200)');
axis tight equal
caxis([0 3])
%%
imagesc(Imap2(:,:,200)');
axis tight equal
caxis([0 4])
%%
Imap2 = zeros(siz2);
Imap2(Imapbla == 1) = 1;
Imap2(Imapbla == 2) = 1;

Imap2(ImapLkid == 1) = 2;
Imap2(ImapLkid == 2) = 2;

Imap2(ImapRkid == 1) = 3;
Imap2(ImapRkid == 2) = 3;

JI2= CalcuJI(Imap2,pM1GT,K-1);
disp(JI2);
%%
%In = pM1E2; InGT = blaGT;
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
%mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = pM1E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);

In = pM1E2; InGT = RkidGT;
mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.01:2.79 2.8];
xlim([0 2.8])

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