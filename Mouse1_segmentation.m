%Tr_M2 M3 M4 Te M1
Xtr = [[M2E2(mask1); M3E2(mask3); M4E2(mask4)] [M2E3(mask1); M3E3(mask3); M4E3(mask4)]...
      [M2E4(mask1); M3E4(mask3); M4E4(mask4)] ];
Xte = [M1E2(mask1) M1E3(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask1); M3GT(mask3); M4GT(mask4)];
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
%初期値
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
%EM1
clearvars atlas
atlas  = atlasfunc1(sig1,sig2,siz,mask1,M2GT,M3GT,M4GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,~]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask1,siz,30);
JI= CalcuJI(Imap,M1GT,K1-1);
disp("EM_MAP result")
disp(JI);

clearvars Xtr XTe XGTtr
%%
%mask
blamask = zeros(siz); blamaxcomp = zeros(siz);
L = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L.PixelIdxList)); 
blamaxcomp(L.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(tmp) = 1;
blamask = logical(and(blamask,mask1)); 

Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L.PixelIdxList)); 
Lmaxcomp(L.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(tmp) = 1;
Lkidmask = logical(and(Lkidmask,mask1));

Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L.PixelIdxList)); 
Rmaxcomp(L.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask(tmp) = 1;
Rkidmask = logical(and(Rkidmask,mask1));
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

clearvars blamaxcomp Rmaxcomp Lmaxcomp L1 XX YY ZZ atlas x y z tmp tmp1 tmp2 tmp3
%%
imagesc(Lkidmask(:,:,200)');
axis tight equal off
%%
Xtebla  = M1E1(blamask);
XteLkid = M1E1(Lkidmask);
XteRkid = M1E1(Rkidmask);
%%
Xtebla  = [M1E3(blamask)  M1E4(blamask) ];
XteLkid = [M1E3(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E3(Rkidmask) M1E4(Rkidmask)];
%%
%Tr_M2 M3 M4 Te M1
Xtebla  = [M1E2(blamask)  M1E3(blamask)  M1E4(blamask)];
XteLkid = [M1E2(Lkidmask) M1E3(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E2(Rkidmask) M1E3(Rkidmask) M1E4(Rkidmask)];
%%
Xtebla  = [M1E1(blamask)  M1E2(blamask)  M1E3(blamask)  M1E4(blamask)];
XteLkid = [M1E1(Lkidmask) M1E2(Lkidmask) M1E3(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E2(Rkidmask) M1E3(Rkidmask) M1E4(Rkidmask)];

%%
%atlas
sig3 = 10;
clearvars atlasbla atlasLkid atlasRkid
atlasbla   = atlasfunc2(sig3,siz,mask1,blamask,GMMpro,0.8,1);
atlasLkid  = atlasfunc2(sig3,siz,mask1,Lkidmask,GMMpro,0.8,2);
atlasRkid  = atlasfunc2(sig3,siz,mask1,Rkidmask,GMMpro,0.8,3);

%初期値
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

%EM2
[Imapbla,~,PPbla,~,~,~,~,~]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,K2,blamask,siz,30);

[ImapLkid,~,PPLkid,~,~,~,~,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,K2,Lkidmask,siz,30);

[ImapRkid,~,PPRkid,~,~,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,K2,Rkidmask,siz,30);

Imap2 = zeros(siz)+4;
Imap2(M1GT==0 ) = 0;
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;

JI = CalcuJI(Imap,M1GT,K1-1);
disp(JI)
JI = CalcuJI(Imap2,M1GT,K1-1);
disp(JI);

clearvars atlasbla atlasLkid atlasRkid Sbla SLkid SRkid Xtebla XteLkid XteRkid
%%
imagesc(M1E2(:,:,200)');
axis tight equal off
colormap(gray)
caxis([0 0.7])
%%
save_raw(Imap2,'C:\\Users\\yourb\\Desktop\\new3\\ImapM1.raw','*uint8');

%%
%mask2
blamask2 = zeros(siz); blamaxcomp = zeros(siz);
L = bwconncomp(Imap2 == 1);
[~,idx] = max(cellfun(@numel,L.PixelIdxList)); 
blamaxcomp(L.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask2(tmp) = 1;
blamask2 = logical(and(blamask2,mask1)); 

Lmaxcomp = zeros(siz); Lkidmask2 = zeros(siz);
L = bwconncomp(Imap2 == 2);
[~,idx] = max(cellfun(@numel,L.PixelIdxList)); 
Lmaxcomp(L.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask2(tmp) = 1;
Lkidmask2 = logical(and(Lkidmask2,mask1));

Rmaxcomp = zeros(siz); Rkidmask2 = zeros(siz);
L = bwconncomp(Imap2 == 3);
[~,idx] = max(cellfun(@numel,L.PixelIdxList)); 
Rmaxcomp(L.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask2(tmp) = 1;
Rkidmask2 = logical(and(Rkidmask2,mask1));
LRAND = and(Rkidmask2,Lkidmask2);
LRkidmask = or(Rkidmask2,Lkidmask2);

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
Rkidmask2 = logical(Rkidmask2 - LRAND + tmp2); 

tmp = zeros(siz);
tmp(RRkid > RLkid) = 1;
tmp2 = and(tmp,LRAND);
Lkidmask2 = logical(Lkidmask2 - LRAND + tmp2);

clearvars blamaxcomp Rmaxcomp Lmaxcomp L1 XX YY ZZ atlas x y z tmp tmp1 tmp2 tmp3
%%
save_raw(Imap2,'C:\\Users\\yourb\\Desktop\\new3\\ImapM1.raw','*uint8');

%%
%Reaginal term bladder
clearvars PPout GraphModel

PPtemp1 = zeros(siz); PPtemp1(blamask) = PPbla(:,1)+PPbla(:,2);
PPtemp1 = imgaussfilt3(PPtemp1,5);
PPtemp2 = 1.0 - PPtemp1;
PPout(:,1) = PPtemp1(blamask2);
PPout(:,2) = PPtemp2(blamask2);
PPout = -log(PPout+eps);

Im = M1E2(blamask2);
GraphModel = CreateFullyConnectedGraphWithMask(blamask2);

clearvars PPtemp1 PPtemp2 PPorgannew
%%
[sigma,lambda] =ndgrid(0.001:0.005:0.1,0.01:0.25:2);
lambda = lambda(:);
sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);
%sigma = 0.005; lambda = 0.01;
%%
%GC bla
clearvars  edgeWeight  Bound
edgeWeight(:,1) = GraphModel.Hi;
edgeWeight(:,2) = GraphModel.Hj;
Z = (Im(GraphModel.Hi)-Im(GraphModel.Hj)).^2;

sigma = 0.006; lambda = 0.26; n =1;
%n= 1;
%for n = 1:160
    terminalWeights = lambda(n) .* PPout;
    Bound = exp(-Z ./ (2*sigma(n)^2)) ./ GraphModel.dist;

    edgeWeight(:,3) = Bound;
    edgeWeight(:,4) = Bound;
    
    [cut, labels] = graphCutMex(terminalWeights,edgeWeight);
    
    Outputtem = zeros(siz);
    Outputtem(blamask2) = labels;
    JI= CalcuJI(Outputtem,M1GT,1);
    disp(JI);
    OutputJI(n) = JI;
    disp(n);
%end
OutGC(Outputtem == 1) = 1;
%%
slice = 70;
subplot(1,3,1);
imagesc(OutGC(:,:,slice)');
axis tight equal off

subplot(1,3,2);
imagesc(Outputtem(:,:,slice)');
axis tight equal off

subplot(1,3,3);
imagesc(M1GT(:,:,slice)');
axis tight equal off
%%
JI = CalcuJI(Imap2,M1GT,3);
disp(JI)













%%
%Reaginal term L.kidney
clearvars PPout GraphModel

PPtemp1 = zeros(siz); PPtemp1(Lkidmask) = PPLkid(:,1)+PPLkid(:,2);
PPtemp1 = imgaussfilt3(PPtemp1,5);
PPtemp2 = 1.0 - PPtemp1;

PPout(:,1) = PPtemp1(Lkidmask2);
PPout(:,2) = PPtemp2(Lkidmask2);
PPout = -log(PPout+eps);

Im = M1E2(Lkidmask2);
GraphModel = CreateFullyConnectedGraphWithMask(Lkidmask2);

clearvars PPtemp1 PPtemp2 PPorgannew
%%
shpepri = Iw(Lkidmask2); 
te = (1-(shpepri(GraphModel.Hi)-shpepri(GraphModel.Hj))./GraphModel.dist)./2;
te =real(sqrt(te));

%%
%GC Lkid
clearvars edgeWeight terminalWeights

edgeWeight(:,1) = GraphModel.Hi;
edgeWeight(:,2) = GraphModel.Hj;
Z = (Im(GraphModel.Hi)-Im(GraphModel.Hj)).^2;
sigma = 0.0019; lambda = 0.039; n =1; c = 1;

%for n = 1:525
    terminalWeights = lambda(n) .* PPout;
    Bound = exp(-Z ./ (2*sigma(n)^2)) ./ GraphModel.dist;
    %edgeWeight(:,3) = c(n)*Bound + (1-c(n))*te;
    edgeWeight(:,3) = Bound;
    edgeWeight(:,4) =   edgeWeight(:,3);
    
    [~, labels] = graphCutMex(terminalWeights,edgeWeight);
    
    Outputtem = zeros(siz);
    Outputtem(Lkidmask2) = labels;
    JI= CalcuJI(Outputtem,M1GT-1,1);
    OutputJI(n) = JI;
    disp(JI); disp(n);
%end

OutGC2(Outputtem == 1) = 2;
%%
slice = 200;
subplot(1,3,1);
imagesc(Imap2(:,:,slice)');
axis tight equal off

subplot(1,3,2);
imagesc(Outputtem(:,:,slice)');
axis tight equal off

subplot(1,3,3);
imagesc(M1GT(:,:,slice)');
axis tight equal off









%%
%Reaginal term R.kidney
clearvars PPout GraphModel

PPtemp1 = zeros(siz); PPtemp1(Rkidmask) = PPRkid(:,1)+PPRkid(:,2);
PPorgannew = zeros(siz);
PPorgannew(Rkidmask2) = PPtemp1(Rkidmask2);
PPtemp1 = imgaussfilt3(PPorgannew,5);
PPtemp2 = 1.0 - PPtemp1;

PPout(:,1) = PPtemp1(Rkidmask2);
PPout(:,2) = PPtemp2(Rkidmask2);
PPout = -log(PPout+eps);
Im = M1E2(Rkidmask2);

GraphModel = CreateFullyConnectedGraphWithMask(Rkidmask2);
clearvars PPtemp1 PPtemp2 PPorgannew
%%
shpepri = Iw(Rkidmask2); 
te = (1-(shpepri(GraphModel.Hi)-shpepri(GraphModel.Hj))./GraphModel.dist)./2;
te =real(sqrt(te));
%%
[sigma,lambda] =ndgrid(0.0001:0.0002:0.002,0.041:0.002:0.05);
lambda = lambda(:);
sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);
%%
%GC Rkid
clearvars edgeWeight terminalWeights

edgeWeight(:,1) = GraphModel.Hi;
edgeWeight(:,2) = GraphModel.Hj;
Z =   (Im(GraphModel.Hi)-Im(GraphModel.Hj)).^2;
sigma = 0.0013; lambda = 0.013; n =1; c = 1.0;

%for n = 1:525
    terminalWeights = lambda(n) .* PPout;
    Bound = exp(-Z ./ (2*sigma(n)^2)) ./ GraphModel.dist;
    %edgeWeight(:,3) = c(n)*Bound + (1-c(n))*te;
    edgeWeight(:,3) = Bound;
    edgeWeight(:,4) =   edgeWeight(:,3);
    
    [~, labels] = graphCutMex(terminalWeights,edgeWeight);
    
    Outputtem = zeros(siz); 
    Outputtem(Rkidmask2) = labels;
    JI= CalcuJI(Outputtem,M1GT-2,1);
    disp(JI);
    OutputJI(n) = JI;
    disp(n);
%end
OutGC2(Outputtem == 1) = 3;

%%
imagesc(Outputtem(:,:,200)');
axis tight equal off
colormap(map)
caxis([0 4])
%%
save_raw(OutGC2,'C:\\Users\\yourb\\Desktop\\new3\\GC2M1.raw','*uint8');