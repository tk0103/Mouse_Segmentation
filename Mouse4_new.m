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
%%
imagesc(Imap2(:,:,200)');
axis tight equal off
%caxis([0 0.7])
%colormap(gray)

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
out = Imap2;
out(Imap2 == 4) = 0;
%imagesc(out(110:290,130:310,200)');
%imagesc(out(220:350,190:320,73)');
imagesc(out(:,:,200)');
axis tight equal off
colormap(map)
caxis([0 4])
%%
imagesc(wM4E1(110:290,130:310,200)');
%imagesc(M4E2(220:350,190:320,73)');
axis tight equal off
colormap(gray)
caxis([0 0.7])
%%
figure;
%In = M4E2; InGT = blaGT;
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
%mu = GMMMubla; sigma = sqrt(GMMSigmabla);

In = M4E2; InGT = LkidGT;
mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);
% 
% In = M4E2; InGT = RkidGT;
% mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
% mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.005:1.2 1.2];
xlim([0 1.2])

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
edge =[0 0:0.005:1.0 1.0];
In = M4E2; InGT = M4GT;
hold on
%histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
%histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);
%%
figure;

edge =[0 0:0.005:1.0 1.0];
In = M4E2; InGT = M4GT;
mu = SS.mu; sigma = sqrt(SS.Sigma);
mu = GMMMu; sigma = sqrt(GMMSigma);

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);

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
plot(edge,y4,'Color',[128 0 128]/255,'LineWidth',2)
hold off



%%
blamask2 = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap2 == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask2(tmp) = 1;
blamask2 = logical(and(blamask2,mask4)); 

Lmaxcomp = zeros(siz); Lkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask2(tmp) = 1;
Lkidmask2 = logical(and(Lkidmask2,mask4));

Rmaxcomp = zeros(siz); Rkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask2(tmp) = 1;
Rkidmask2 = logical(and(Rkidmask2,mask4));
LRAND = and(Rkidmask2,Lkidmask2);

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

Output = zeros(siz);

%%
%GC bladder
clearvars PPout GraphModel
masktm = blamask;
masktm2 = blamask2;
GT = zeros(siz); GT(M4GT == 1) = 1;
PPorgan = zeros(siz); PPorgan(masktm) = PPbla(:,1) + PPbla(:,2);
PPorgan = imgaussfilt3(PPorgan,5);
PPback = 1.0 - PPorgan;
PPback(PPback <0 ) = 0;
PPout(:,1) = PPorgan(masktm2); PPout(:,2) = PPback(masktm2);

PPout = -log(PPout+eps);
GraphModel = CreateFullyConnectedGraphWithMask(masktm2);

Kmat = [1,1; 1,1]; K3 = 2;
sumIm = (M4E1+M4E2)/2;
sumIm = sumIm(masktm2);

%%
%Reaginal term
clearvars PPout GraphModel
masktm = Lkidmask;
masktm2 = Lkidmask2;
GT = zeros(siz); GT(M4GT == 2) = 1;
PPorgan = zeros(siz); PPorgan(masktm) = (PPLkid(:,1)+PPLkid(:,2));
PPorgannew = zeros(siz);
PPorgannew(masktm2) = PPorgan(masktm2);
PPorgan = imgaussfilt3(PPorgannew,5);
PPback = 1.0 - PPorgan;
PPout(:,1) = PPorgan(masktm2); PPout(:,2) = PPback(masktm2);
PPout = -log(PPout+eps);

Kmat = [1,1; 1,1]; K3 = 2;
sumIm = (M4E1 + M4E2 )/2;
sumIm = sumIm(masktm2);

GraphModel = CreateFullyConnectedGraphWithMask(masktm2);
%%
%Reaginal term
clearvars PPout GraphModel
masktm = Rkidmask;
masktm2 = Rkidmask2;
GT = zeros(siz); GT(M4GT == 3) = 1;
PPorgan = zeros(siz); PPorgan(masktm) = PPRkid(:,1)+PPRkid(:,2);
PPorgannew = zeros(siz);
PPorgannew(masktm2) = PPorgan(masktm2);
PPorgan = imgaussfilt3(PPorgannew,5);
PPback = 1.0 - PPorgan;
PPout(:,1) = PPorgan(masktm2); PPout(:,2) = PPback(masktm2);
PPout = -log(PPout+eps);

Kmat = [1,1; 1,1]; K3 = 2;
sumIm = (M4E1 + M4E2 )/2;
sumIm = sumIm(masktm2);

GraphModel = CreateFullyConnectedGraphWithMask(masktm2);
%%
[sig,lambda] =ndgrid(0.00005:0.005:0.2,0.2:0.1:0.4);
lambda = lambda(:);
sig = sig(:);
OutputJI = zeros(size(sig,1),1);
%%
%�t��
%for n = 1:size(sig,1)
n = 1;
lambda = 0.5;
sig = 0.00505;
N = size(PPout,1);
CurLabel = zeros(N,1)+2;
PreLabel = zeros(N,1);
Outtemp = zeros(siz);
PropLabel = ones(N,1);
flag = 0;

while(flag~=1)
    
    GraphModel.Vs = PPout((1:N)'+(CurLabel(:)-1)*N)*lambda(n);
    GraphModel.Vt = PPout((1:N)'+(PropLabel(:)-1)*N)*lambda(n);
    GraphModel = SetNWeight_binary(GraphModel, sumIm, CurLabel,PropLabel,sig(n), Kmat);
    
    [lowerBound, label] = qpboMex([GraphModel.Vs,GraphModel.Vt],...
        [GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
    
    label = logical(label);
    CurLabel(label) = PropLabel(label);
        
    GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
    Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H01);
    disp(E);
    
    
    if CurLabel == PreLabel
        flag = 1;
        Outtemp(masktm2) = CurLabel;
    end
    
    PreLabel = CurLabel;
    
    PreE = E;
end

JI3= CalcuJI(Outtemp,GT,1);
disp(JI3);
OutputJI(n) = JI3';
disp(n);
%end
Output(Outtemp == 1) = 3;
%%
save_raw(Outtemp,'C:\\Users\\yourb\\Desktop\\M3GC.raw','*uint8')
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

%%
aaa = zeros(siz); aaaa(Rkidmask2) = label;
%%
imagesc(Output(:,:,73)');
axis tight equal off
caxis([0 4])
colormap(map)
%%
JI2= CalcuJI(Output,M4GT,K2);
disp(JI2);
%%
clearvars GraphModel
GraphModel = CreateFullyConnectedGraphWithMask(masktm2);
   %%
%for n = 1:size(sig,1)
n =1;
lambda = 5;
sig = 0.05;
    N = size(PPout,1);
    CurLabel = zeros(N,1)+2;
    PreLabel = zeros(N,1);
    Outtemp = zeros(siz);
    PropLabel = ones(N,1);
    
    flag = 0; PreE = 0; 
 
    while(flag~=1)
        GraphModel.Vs = PPout((1:N)'+(CurLabel(:)-1)*N)*lambda(n);
        GraphModel.Vt = PPout((1:N)'+(PropLabel(:)-1)*N)*lambda(n);
        GraphModel = SetNWeight_binary(GraphModel, sumIm, CurLabel,PropLabel,sig(n), Kmat);
  
        [lowerBound, label] = qpboMex([GraphModel.Vs,GraphModel.Vt],...
            [GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]); 
      
        label = logical(label);
        CurLabel(label) = PropLabel(label);
        
        GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
        Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H01);
        E = Eunary + Epairwise;
        disp(E);
       
        if CurLabel == PreLabel
            flag = 1;
            Outtemp(masktm2) = CurLabel;
        end
        
        PreLabel = CurLabel;
        
        PreE = E;
    end
    
    JI3= CalcuJI(Outtemp,GT,1);
    disp(JI3);
    OutputJI(n) = JI3';
    disp(n);
%end

Output(Outtemp == 1) = 1;

%%
%Reaginal term
tempPP1 = zeros(siz); tempPP2 = zeros(siz); tempPP3 = zeros(siz);
tempPP1(pmask1) = PP(:,1); tempPP2(pmask1) = PP(:,2); tempPP3(pmask1) = PP(:,3);
tempPP1 = imgaussfilt3(tempPP1,5);
tempPP2 = imgaussfilt3(tempPP2,5);
tempPP3 = imgaussfilt3(tempPP3,5);
tempPP4 = 1 - tempPP1 - tempPP2 - tempPP3;
PPout(:,1) = tempPP1(pmask1); PPout(:,2) = tempPP2(pmask1); 
PPout(:,3) = tempPP3(pmask1); PPout(:,4) = tempPP4(pmask1);

RP = cell(1,K);
for k = 1:K
    RP{1,k} = -log(PP(:,k)+eps);
    RP{1,k}(isnan(RP{1,k})) = 0;
end
RP = cell2mat(RP);
RP = real(RP);

%GraphCut
GraphModel = CreateFullyConnectedGraphWithMask(pmask1);



%Shape term
[E1,E2] = Create_shape(L,GraphModel,pmask1,siz2);

%voronoi
voronoiIn = zeros(siz2);
for k = 1:K-1
    LL = bwconncomp(Imap==k);
    numPixels = cellfun(@numel,LL.PixelIdxList);
    [~,idx] = max(numPixels);
    voronoiIn(LL.PixelIdxList{idx}) = k;
end
voronoiFig = zeros(siz2);
[voronoiOut,~] = mistVoronoiDistanceTransform(uint8(voronoiIn(pmask1)));
voronoiFig(pmask1) = voronoiOut;

%%
% for n = 1:27
n = 1;
N = size(RP,1);
CurLabel = zeros(N,1)+K;
PreLabel = zeros(N,1);
Output = zeros(siz2);
flag = 0;
PreE = 0;

%%
Sigmat =  abs(bsxfun(@minus,GMMMu(:,1),GMMMu(:,1)'))*h(n) + eye(K);
PropLabel = double(voronoiOut);
PropLabel(PropLabel == 0) = 1;

while(flag ~=1)
    GraphModel = SetTWeights(GraphModel,RP,CurLabel,PropLabel,lambda(n),N);
    GraphModel = SetNWeights(GraphModel,pM1E2(pmask1),CurLabel,PropLabel,Sigmat,graydiff,shape,E1,E2,c(n));
    [lowerBound, labels] = qpboMex([GraphModel.Vs,GraphModel.Vt],[GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
    labels = logical(labels);
    CurLabel(labels) = PropLabel(labels);
    
    GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
    Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H00);
    E = Eunary + Epairwise;
    
    disp(E);
    if CurLabel == PreLabel
        flag = 1;
        Output(pmask1) = CurLabel;
 
    end
    
    PreLabel = CurLabel;
    PreE = E;
end

JI= CalcuJI(Output,pM1GT,K-1);
disp("GraphCut_JI")
disp(JI);
% 
OutputJI(n,:) = JI';
disp(n);
% end
%%
save_raw(Output,'C:\\Users\\yourb\\Desktop\\M1GC.raw','*uint8')
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
%%
Imapout = Imap;
Imapout(Imap==4) = 0;
pM1GTout = pM1GT;
pM1GTout(pM1GT==4) = 0;
%%
imagesc(Imap2(:,:,200)');
%%
slice1 = 206;
slice2 = 66;

subplot(2,2,1)
imagesc(Imapout(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,2)
imagesc(pM1GTout(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,3)
imagesc(Imapout(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,4)
imagesc(pM1GTout(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

%%
blamask2 = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap2 == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask2(tmp) = 1;
blamask2 = logical(and(blamask2,mask4)); 

Lmaxcomp = zeros(siz); Lkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask2(tmp) = 1;
Lkidmask2 = logical(and(Lkidmask2,mask4));

Rmaxcomp = zeros(siz); Rkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask2(tmp) = 1;
Rkidmask2 = logical(and(Rkidmask2,mask4));
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
%%
imagesc(voronoi(:,:,200)');
axis tight equal off

%%
%Reaginal term
clearvars PPout
masktm = LRkidmask;

PPtemp1 = zeros(siz); PPtemp1(Lkidmask) = PPLkid(:,1)+PPLkid(:,2);
PPorgannew = zeros(siz);
PPorgannew(masktm) = PPtemp1(masktm);
PPtemp1 = imgaussfilt3(PPorgannew,5);

PPtemp2 = zeros(siz); PPtemp2(Rkidmask) = PPRkid(:,1)+PPRkid(:,2);
PPorgannew = zeros(siz);
PPorgannew(masktm) = PPtemp2(masktm);
PPtemp2 = imgaussfilt3(PPorgannew,5);
PPtemp3 = 1.0 - PPtemp1 - PPtemp2;

PPout(:,1) = PPtemp1(masktm);
PPout(:,2) = PPtemp2(masktm);
PPout(:,3) = PPtemp3(masktm);
PPout = -log(PPout+eps);

sumIm = M4E2;
sumIm = sumIm(masktm);

GraphModel = CreateFullyConnectedGraphWithMask(masktm);

clearvars PPtemp1 PPtemp2 PPtemp3 PPorgannew

%%
%voronoi
voronoi = zeros(siz);
voronoi(Lkidmask2) = 1;
voronoi(Rkidmask2) = 2;
%%
lambda = 1;
h = 1;
n = 1; 

flag = 0; PreE = 0;
N = size(PPout,1);
CurLabel = zeros(N,1)+K1-1;
PreLabel = zeros(N,1);
Output = zeros(siz);
PropLabel = double(voronoi(LRkidmask));
PropLabel(PropLabel == 0) = 1;

newGMMMMu(1,:) = GMMMuLkid(1,:);
newGMMMMu(2,:) = GMMMuRkid(1,:);
newGMMMMu(3,:) = (GMMMuLkid(3,:) + GMMMuRkid(3,:))./2;
Sigmat =  abs(bsxfun(@minus,newGMMMMu(:,1),newGMMMMu(:,1)'))*h(n) + eye(K1-1);
Kmat = ones(K1-1);

while(flag ~=1)
    GraphModel = SetTWeights(GraphModel,PPout,CurLabel,PropLabel,lambda(n),N);
    GraphModel = SetNWeights(GraphModel,sumIm,CurLabel,PropLabel,Sigmat,Kmat);
    [lowerBound, labels] = qpboMex([GraphModel.Vs,GraphModel.Vt],[GraphModel.Hi,...
        GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
    labels = logical(labels);
    CurLabel(labels) = PropLabel(labels);
    
    GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
    Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H00);
    E = Eunary + Epairwise;
    disp(E);
    
    if CurLabel == PreLabel
        flag = 1;
        Output(masktm) = CurLabel;
    end
    
    PreLabel = CurLabel;
    PreE = E;
end

JI= CalcuJI(Output,M4GT-1,K1-1);
disp("GraphCut_JI")
disp(JI);

%%
subplot(1,2,1);
imagesc(Imap2(:,:,200)');
axis tight equal 
caxis([0 4])

subplot(1,2,2);
imagesc(Output(:,:,200)');
axis tight equal 
caxis([0 4])

%%
OutputJI(n,:) = JI';
disp(n);
% end
%%
[sig,lambda] =ndgrid(0.00005:0.005:0.2,0.2:0.1:0.4);
lambda = lambda(:);
sig = sig(:);
OutputJI = zeros(size(sig,1),1);


imagesc(Output(:,:,200)');
