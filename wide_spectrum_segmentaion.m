mask = mask4;
GT  =M4GT;
%%
%train_mouse2_mouse3 mouse4 test_mouse1
Xtr = [wM2E1(mask2); wM3E1(mask3); wM4E1(mask4) ];
Xte = wM1E1(mask1);
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];

%%
%train_mouse1 mouse3 mouse4 test_mouse2
Xtr = [wM1E1(mask1); wM3E1(mask3); wM4E1(mask4) ];
Xte = wM2E1(mask2);
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];

%%
%train_mouse1 mouse2 mouse4 test_mouse3
Xtr = [wM1E1(mask1); wM2E1(mask2); wM4E1(mask4) ];
Xte = wM3E1(mask3);
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];

%%
%train_mouse1 mouse2 mouse3 test_mouse4
Xtr = [wM1E1(mask1); wM2E1(mask2); wM3E1(mask3) ];
Xte = wM4E1(mask4);
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask4)];
%%
%initial_value
clearvars SS
for k = 1:K1
    tmp1 = Xtr(:,1);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.Sigma(:,:,k) = cov(tmp1(XGTtr == k));
end

%%
%Atlas_guided EM,Mouse1
clearvars atlas
atlas  = atlasfunc1(sig1,sig2,K1,siz,mask1,M2GT,M3GT,M4GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,~]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask1,siz,30);
JI= CalcuJI(Imap,M1GT,K1-1);
disp("EM_MAP result")
disp(JI);
%%
%Atlas_guided EM,Mouse2
clearvars atlas
atlas  = atlasfunc1(sig1,sig2,K1,siz,mask2,M1GT,M3GT,M4GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,~] = ...
    AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask2,siz,30);
JI= CalcuJI(Imap,M2GT,K1-1);
disp("EM_MAP result")
disp(JI);
%%
%Atlas_guided EM,Mouse3
clearvars atlas
atlas  = atlasfunc1(sig1,sig2,K1,siz,mask3,M1GT,M2GT,M4GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,~] ...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask3,siz,30);
JI= CalcuJI(Imap,M3GT,K1-1);
disp("EM_MAP result")
disp(JI);

%%
%EM Mosue4
clearvars atlas
atlas  = atlasfunc1(sig1,sig2,K1,siz,mask4,M1GT,M2GT,M3GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,~]...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask4,siz,30);
JI= CalcuJI(Imap,M4GT,K1-1);
disp("EM_MAP result")
disp(JI);

%%
blamask = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(tmp) = 1;
blamask = logical(and(blamask,mask1)); 

Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L1 = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(tmp) = 1;
Lkidmask = logical(and(Lkidmask,mask1));

Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L1 = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
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
%Mouse1 Test
Xtebla  = wM1E1(blamask);
XteLkid = wM1E1(Lkidmask);
XteRkid = wM1E1(Rkidmask);

%%
%Mouse2 Test
Xtebla  = wM2E1(blamask);
XteLkid = wM2E1(Lkidmask);
XteRkid = wM2E1(Rkidmask);

%%
%Mouse3 Test
Xtebla  = wM3E1(blamask);
XteLkid = wM3E1(Lkidmask);
XteRkid = wM3E1(Rkidmask);

%%
%Mouse4 Test
Xtebla  = wM4E1(blamask);
XteLkid = wM4E1(Lkidmask);
XteRkid = wM4E1(Rkidmask);
%%
sig3 = 10;
clearvars atlasbla atlasLkid atlasRkid
atlasbla   = atlasfunc2(sig3,siz,mask,blamask,GMMpro,0.8,1);
atlasLkid  = atlasfunc2(sig3,siz,mask,Lkidmask,GMMpro,0.2,2);
atlasRkid  = atlasfunc2(sig3,siz,mask,Rkidmask,GMMpro,0.2,3);

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

[ImapLkid,~,PPLkid,GMMMuLkid,GMMSigmaLkid,~,FeatL,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,K2,Lkidmask,siz,30);

[ImapRkid,~,PPRkid,GMMMuRkid,GMMSigmaRkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,K2,Rkidmask,siz,30);

Imap2 = zeros(siz)+4;
Imap2(GT==0 ) = 0;
Imap2(Imapbla == 1) = 1; Imap2(Imapbla == 2) = 1;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;

clearvars atlasbla atlasLkid atlasRkid Sbla SLkid SRkid Xtebla XteLkid XteRkid

JI = CalcuJI(Imap,GT,K1-1);
disp(JI);
JI = CalcuJI(Imap2,GT,K1-1);
disp(JI);
%%
save_raw(Imap2,'C:\\Users\\yourb\\Desktop\\new3\\ImapWM4.raw','*uint8');
%%
blamask2 = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap2 == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask2(tmp) = 1;
blamask2 = logical(and(blamask2,mask)); 

Lmaxcomp = zeros(siz); Lkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask2(tmp) = 1;
Lkidmask2 = logical(and(Lkidmask2,mask));

Rmaxcomp = zeros(siz); Rkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask2(tmp) = 1;
Rkidmask2 = logical(and(Rkidmask2,mask));
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

clearvars blamaxcomp Rmaxcomp Lmaxcomp L1 XX YY ZZ atlas x y z tmp tmp1 tmp2 tmp3 Xtr XGTtr Xte
%%
imagesc(OutGC(:,:,200)');
axis tight equal off
colormap(map)
caxis([0 4])
%%
save_raw(Imap2,'C:\\Users\\yourb\\Desktop\\new3\\Imap2_Mouse2_wide.raw','*uint8');
%%
%Reaginal term bladder
clearvars PPout GraphModel

PPtemp1 = zeros(siz); PPtemp1(blamask) = PPbla(:,1)+PPbla(:,2);
PPorgannew = zeros(siz);
PPorgannew(blamask2) = PPtemp1(blamask2);
PPtemp1 = imgaussfilt3(PPorgannew,5);
PPtemp2 = 1.0 - PPtemp1;

PPout(:,1) = PPtemp1(blamask2);
PPout(:,2) = PPtemp2(blamask2);
PPout = -log(PPout+eps);

Im = wM4E1(blamask2);
GraphModel = CreateFullyConnectedGraphWithMask(blamask2);

clearvars PPtemp1 PPtemp2 PPorgannew

%δNγχ
[sigma,lambda] =ndgrid(0.001:0.005:0.1,0.01:0.25:2);
lambda = lambda(:);
sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);
%sigma = 0.005; lambda = 0.01;
%%
%bladder
clearvars edgeWeight terminalWeights
edgeWeight(:,1) = GraphModel.Hi;
edgeWeight(:,2) = GraphModel.Hj;
Z = (Im(GraphModel.Hi)-Im(GraphModel.Hj)).^2;
%sigma = 0.011; lambda = 3.01; n =1;

for n = 1:160
    terminalWeights = lambda(n) .* PPout;
    Bound = exp(- Z ./ (2*sigma(n)^2)) ./ GraphModel.dist;
    edgeWeight(:,3) = Bound;
    edgeWeight(:,4) = Bound;
    
    [~, labels] = graphCutMex(terminalWeights,edgeWeight);
    
    Outputtem = zeros(siz);
    Outputtem(blamask2) = labels;
    JI= CalcuJI(Outputtem,GT,1);
    disp(JI);
    OutputJI(n) = JI;
    disp(n);
end
OutGC(Outputtem == 1) = 1;
%%
slice = 80;
subplot(1,3,1);
imagesc(Imap2(:,:,slice)');
axis tight equal off

subplot(1,3,2);
imagesc(OutGC(:,:,slice)');
axis tight equal off

subplot(1,3,3);
imagesc(GT(:,:,slice)');
axis tight equal off
%%
JI = CalcuJI(Imap2,M2GT,3);
disp(JI)


%%
%Reaginal term L.kidney
clearvars PPout GraphModel

PPtemp1 = zeros(siz); PPtemp1(Lkidmask) = PPLkid(:,1)+PPLkid(:,2);
PPorgannew = zeros(siz);
PPorgannew(Lkidmask2) = PPtemp1(Lkidmask2);
PPtemp1 = imgaussfilt3(PPorgannew,5);
PPtemp2 = 1.0 - PPtemp1;

PPout(:,1) = PPtemp1(Lkidmask2);
PPout(:,2) = PPtemp2(Lkidmask2);
PPout = -log(PPout+eps);

Im = wM4E1(Lkidmask2);
GraphModel = CreateFullyConnectedGraphWithMask(Lkidmask2);

clearvars PPtemp1 PPtemp2 PPorgannew
%%
[sigma,lambda] =ndgrid(0.0001:0.0002:0.002,0.001:0.002:0.04);
lambda = lambda(:);
sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);
%%
clearvars edgeWeight terminalWeights
edgeWeight(:,1) = GraphModel.Hi;
edgeWeight(:,2) = GraphModel.Hj;
Z = (Im(GraphModel.Hi)-Im(GraphModel.Hj)).^2;
%sigma = 0.0017; lambda = 0.019; n =1;

for n = 1:200
    terminalWeights = lambda(n) .* PPout;
    Bound = exp(- Z ./ (2*sigma(n)^2)) ./ GraphModel.dist;
    edgeWeight(:,3) = Bound;
    edgeWeight(:,4) = Bound;
    
    [~, labels] = graphCutMex(terminalWeights,edgeWeight);
    
    Outputtem = zeros(siz);
    Outputtem(Lkidmask2) = labels;
    JI= CalcuJI(Outputtem,GT-1,1);
    disp(JI);
    OutputJI(n) = JI;
    disp(n);    
end
OutGC(Outputtem ==1) = 2;
%%
imagesc(Outputtem(:,:,200)');

%%
Imap2(Imap2 == 4) =0;
slice = 200;
imagesc(Imap2(:,:,slice)');
axis tight equal off
colormap(map)
caxis([ 0 4])
%%
slice = 200;
subplot(1,3,1);
imagesc(Imap2(:,:,slice)');
axis tight equal off
%%

subplot(1,3,2);
imagesc(Outputtem(:,:,slice)');
axis tight equal off

subplot(1,3,3);
imagesc(GT(:,:,slice)');
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

Im = wM1E1(Rkidmask2);
GraphModel = CreateFullyConnectedGraphWithMask(Rkidmask2);
clearvars PPtemp1 PPtemp2 PPorgannew
%%
[sigma,lambda] =ndgrid(0.0001:0.0002:0.002,0.0001:0.0005:0.001);
lambda = lambda(:);
sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);
%%
clearvars edgeWeight terminalWeights
edgeWeight(:,1) = GraphModel.Hi;
edgeWeight(:,2) = GraphModel.Hj;
Z =   (Im(GraphModel.Hi)-Im(GraphModel.Hj)).^2;
%sigma = 0.0001; lambda = 0.031; n =1;

for n = 1:20
    terminalWeights = lambda(n) .* PPout;
    Bound = exp(-Z ./ (2*sigma(n)^2)) ./ GraphModel.dist;
    edgeWeight(:,3) = Bound;
    edgeWeight(:,4) = Bound;
    
    [~, labels] = graphCutMex(terminalWeights,edgeWeight);
    
    Outputtem = zeros(siz);
    Outputtem(Rkidmask2) = labels;
    JI= CalcuJI(Outputtem,GT-2,1);
    disp(JI);
    OutputJI(n) = JI;
    disp(n);    
end

OutGC(Outputtem ==1) = 3;
%%
imagesc(Imap2(:,:,70)');
axis tight equal off
colormap(map)
caxis([0 4])
%%
save_raw(OutGC,'C:\\Users\\yourb\\Desktop\\new3\\GC_Mouse4_wide.raw','*uint8');