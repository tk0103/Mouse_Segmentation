%train_mouse2_mouse3 test_mouse1
Xtr = [[M1E2(mask1); M2E2(mask2); M3E2(mask3)] [M1E3(mask1); M2E3(mask2); M3E3(mask3)]...
      [M1E4(mask1); M2E4(mask2); M3E4(mask3)] ];
Xte = [M4E2(mask4) M4E3(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
XGTte = M4GT(mask4);
%%
clearvars SS
for k = 1:K1 
    tmp1 = Xtr(:,1); 
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.Sigma(:,:,k) = cov((tmp1(XGTtr == k)));
end 
%%
clearvars SS
for k = 1:K1 
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k)]));
end 
%%
clearvars SS
for k = 1:K1 
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end 

%%
clearvars SS
for k = 1:K1 
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3); tmp4 = Xtr(:,4);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.mu(k,4) = mean(tmp4(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k),tmp4(XGTtr == k)]));
end 
%%
clearvars atlas
atlas  = atlasfunc1(sig1,sig2,K1,siz,mask4,M1GT,M2GT,M3GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,likelihood]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask4,siz,30);
JI1= CalcuJI(Imap,M4GT,K1-1);
disp("EM_MAP result")
disp(JI1);
%%
imagesc(M3E1(:,:,70)');
axis tight equal
%%
imagesc(Imap(:,:,200)');
axis tight equal
caxis([0 4])
%%
blamask = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(tmp) = 1;
blamask = logical(and(blamask,mask4)); 

Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L1 = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L1 = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;

tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(tmp) = 1;
Lkidmask = logical(and(Lkidmask,mask4));

tmp = bwdist(logical(Rmaxcomp)) < power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask(tmp) = 1;
Rkidmask = logical(and(Rkidmask,mask4));
LRAND = and(Rkidmask,Lkidmask);

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

tmp = zeros(siz);
tmp(RRkid > RLkid) = 1;
tmp2 = and(tmp,LRAND);
Lkidmask = logical(Lkidmask - LRAND + tmp2);

%%
imagesc(Rkidmask(:,:,230)');
axis tight equal off
%%
Xtebla  = M4E4(blamask);
XteLkid = M4E4(Lkidmask);
XteRkid = M4E4(Rkidmask);
%%
Xtebla  = [M4E3(blamask)  M4E4(blamask) ];
XteLkid = [M4E3(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E3(Rkidmask) M4E4(Rkidmask)];
%%
Xtebla  = [M4E2(blamask)  M4E3(blamask)  M4E4(blamask)];
XteLkid = [M4E2(Lkidmask) M4E3(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E2(Rkidmask) M4E3(Rkidmask) M4E4(Rkidmask)];
%%
Xtebla  = [M4E1(blamask)  M4E2(blamask)  M4E3(blamask)  M4E4(blamask)];
XteLkid = [M4E1(Lkidmask) M4E2(Lkidmask) M4E3(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E2(Rkidmask) M4E3(Rkidmask) M4E4(Rkidmask)];

%%
GT = zeros(siz);  GT(blamask) = M4GT(blamask);
blaGT = zeros(siz); blaGT(blamask) = 3;
blaGT(GT == 1) = 1;

GT = zeros(siz);  GT(Lkidmask) = M4GT(Lkidmask);
LkidGT = zeros(siz); LkidGT(Lkidmask) = 3;
LkidGT(GT == 2) = 1;

GT = zeros(siz);  GT(Rkidmask) = M4GT(Rkidmask);
RkidGT = zeros(siz); RkidGT(Rkidmask) = 3;
RkidGT(GT == 3) = 1;
%%
GT = zeros(siz);  GT(blamask) = cutM4GT(blamask);
blaGT = zeros(siz); blaGT(blamask) = 3;
blaGT(GT == 1) = 1; blaGT(GT == 5) = 2; 

GT = zeros(siz);  GT(Lkidmask) = cutM4GT(Lkidmask);
LkidGT = zeros(siz); LkidGT(Lkidmask) = 3;
LkidGT(GT == 2) = 1; LkidGT(GT == 6) = 2;

GT = zeros(siz);  GT(Rkidmask) = cutM4GT(Rkidmask);
RkidGT = zeros(siz); RkidGT(Rkidmask) = 3;
RkidGT(GT == 3) = 1; RkidGT(GT == 7) = 2; 

%%
tmp1 = M4E2; tmp2 = M4E3; tmp3 = M4E4; mask = blaGT;
for k = 1:K2
    Sbla.mu(k,1) = mean(tmp1(mask == k));
    Sbla.mu(k,2) = mean(tmp2(mask == k));
    Sbla.mu(k,3) = mean(tmp3(mask == k));
    Sbla.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

mask = LkidGT;
for k = 1:K2
    SLkid.mu(k,1) = mean(tmp1(mask == k));
    SLkid.mu(k,2) = mean(tmp2(mask == k));
    SLkid.mu(k,3) = mean(tmp3(mask == k));
    SLkid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end

mask = RkidGT;
for k = 1:K2
    SRkid.mu(k,1) = mean(tmp1(mask == k));
    SRkid.mu(k,2) = mean(tmp2(mask == k));
    SRkid.mu(k,3) = mean(tmp3(mask == k));
    SRkid.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3
%%

clearvars atlasbla atlasLkid atlasRkid
sig3  = 10;
atlasbla   = atlasfunc2(sig3,siz,mask4,blamask,GMMpro,0.8,1);
atlasLkid  = atlasfunc2(sig3,siz,mask4,Lkidmask,GMMpro,0.2,2);
atlasRkid  = atlasfunc2(sig3,siz,mask4,Rkidmask,GMMpro,0.2,3);


clearvars Sbla SLkid SRkid
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

[Imapbla,~,PPbla,GMMMubla,GMMSigmabla,GMMprobla,Featbla,~]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,K2,blamask,siz,30);

[ImapLkid,~,PPLkid,GMMMuLkid,GMMSigmaLkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,K2,Lkidmask,siz,30);

[ImapRkid,~,PPRkid,GMMMuRkid,GMMSigmaRkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,K2,Rkidmask,siz,10);

Imap2 = zeros(siz) + 4;
Imap2(M4GT == 0) = 0;
Imap2(Imapbla == 1) = 1;  Imap2(Imapbla == 2) = 1;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;

JI1= CalcuJI(Imap,M4GT,K2);
disp(JI1);
JI2= CalcuJI(Imap2,M4GT,K2);
disp(JI2);

