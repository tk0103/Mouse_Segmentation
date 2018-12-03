Xtr = [[pM1E2(pmask1); pM2E2(pmask2)] [pM1E3(pmask1); pM2E3(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)] ];
Xte = [pM3E2(pmask3) pM3E3(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];
XGTte = pM3GT(pmask3);
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
atlas  = atlasfunc2(sig1,sig2,K,siz,pmask3,pM1GT,pM2GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,lilelihood] ...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask3,siz,30);
JI1= CalcuJI(Imap,pM3GT,K-1);
disp("EM_MAP result")
disp(JI1);
%%
imagesc(Imap(:,:,155)');
axis tight equal off
caxis([0 4])
%%
blamask = zeros(siz);
temp = Imap == 1; temp2 = zeros(siz);
L1 = bwconncomp(temp); numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels); 
temp2(L1.PixelIdxList{idx}) = 1; 
temp = bwdist(logical(temp2)) <  power(bwarea(temp(:))/4/pi*3,1/3);
blamask(temp) = 1;
blamask = and(blamask,pmask3); blamask = logical(blamask);
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
Lkidmask = and(Lkidmask,pmask3);
Lkidmask = logical(Lkidmask);

Rkidmask = zeros(siz);
temp = bwdist(logical(temp3)) < Lkidregion;
Rkidmask(temp) = 1;
Rkidmask = and(Rkidmask,pmask3);
Rkidmask = logical(Rkidmask);
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
atlasbla   = atlasfunc3(sig1,siz,pmask3,blamask,GMMpro,blaratio,1);
atlasLkid  = atlasfunc3(sig1,siz,pmask3,Lkidmask,GMMpro,Lkidration,2);
atlasRkid  = atlasfunc3(sig1,siz,pmask3,Rkidmask,GMMpro,Rkidration,3);

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
[Imapbla,~,~,GMMMubla,GMMSigmabla,~,~,~]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,3,blamask,siz,30);

[ImapLkid,~,~,GMMMuLkid,GMMSigmaLkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,3,Lkidmask,siz,30);

[ImapRkid,~,~,GMMMuRkid,GMMSigmaRkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,3,Rkidmask,siz,30);
%%
Imap2 = zeros(siz);
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;


Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
JI2= CalcuJI(Imap2,pM3GT,K-1);

disp(JI2);

%%
temp = zeros(siz);
temp(Rkidmask) = pM2E2(Rkidmask);
%%
imagesc(pM3GT(:,:,200)');
axis tight equal off
caxis([0 4])
%%
imagesc(pM3E2(:,:,80)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
In = pM3E2; InGT = blaGT;
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = pM3E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);


%In = pM3E2; InGT = RkidGT;
%mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
%mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

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
save_raw(temp2,'C:\\Users\\yourb\\Desktop\\temp.raw','*uint8');