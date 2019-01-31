%train_mouse1 mouse2 mouse4 test_mouse3
Xtr = [[M1E2(mask1); M2E2(mask2); M4E2(mask4)] [M1E3(mask1); M2E3(mask2); M4E3(mask4)]...
      [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E2(mask3) M3E3(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
XGTte = M3GT(mask3);
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
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k)]));
end
%%
clearvars SS
for k = 1:K1
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); 
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
atlas  = atlasfunc1(sig1,sig2,K1,siz,mask3,M1GT,M2GT,M4GT);

[Imap,~,~,GMMMu,GMMSigma,GMMpro,~,lilelihood] ...
    = AtlasGuidedEM_kubo(Xte,atlas,SS,K1,mask3,siz,30);
JI1= CalcuJI(Imap,M3GT,K1-1);
disp("EM_MAP result")
disp(JI1);
%%
imagesc(Imap(:,:,200)');
axis tight equal off
caxis([0 4])
%%
blamask = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask(tmp) = 1;
blamask = logical(and(blamask,mask3)); 

Lmaxcomp = zeros(siz); Lkidmask = zeros(siz);
L1 = bwconncomp(Imap == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask(tmp) = 1;
Lkidmask = logical(and(Lkidmask,mask3));


Rmaxcomp = zeros(siz); Rkidmask = zeros(siz);
L1 = bwconncomp(Imap == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask(tmp) = 1;
Rkidmask = logical(and(Rkidmask,mask3));
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
Xtebla  = M3E4(blamask);
XteLkid = M3E4(Lkidmask);
XteRkid = M3E4(Rkidmask);
%%
Xtebla  = [M3E3(blamask)  M3E4(blamask) ];
XteLkid = [M3E3(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E3(Rkidmask) M3E4(Rkidmask)];
%%
Xtebla  = [M3E2(blamask)  M3E3(blamask)  M3E4(blamask)  ];
XteLkid = [M3E2(Lkidmask) M3E3(Lkidmask) M3E4(Lkidmask) ];
XteRkid = [M3E2(Rkidmask) M3E3(Rkidmask) M3E4(Rkidmask) ];
%%
Xtebla  = [M3E1(blamask)  M3E2(blamask)  M3E3(blamask)  M3E4(blamask)];
XteLkid = [M3E1(Lkidmask) M3E2(Lkidmask) M3E3(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E2(Rkidmask) M3E3(Rkidmask) M3E4(Rkidmask)];
%%
clearvars atlasbla atlasLkid atlasRkid
sig3 = 10;
atlasbla   = atlasfunc2(sig3,siz,mask3,blamask,GMMpro,0.8,1);
atlasLkid  = atlasfunc2(sig3,siz,mask3,Lkidmask,GMMpro,0.2,2);
atlasRkid  = atlasfunc2(sig3,siz,mask3,Rkidmask,GMMpro,0.2,3);
%%
tmp  =zeros(siz);
tmp(Lkidmask) = atlasLkid(:,1);
%%
imagesc(tmp(:,:,200)');
axis tight equal
caxis([0 1]);
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
tmp1 = M3E2; tmp2 = M3E3; tmp3 = M3E4; mask = blaGT;
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
%%
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

[Imapbla,~,PPbla,GMMMubla,GMMSigmabla,GMMprobla,Featbla,likelihoodbla]...
    = AtlasGuidedEM_kubo(Xtebla,atlasbla,Sbla,K2,blamask,siz,30);

[ImapLkid,~,PPLkid,GMMMuLkid,GMMSigmaLkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteLkid,atlasLkid,SLkid,K2,Lkidmask,siz,30);

[ImapRkid,~,PPRkid,GMMMuRkid,GMMSigmaRkid,~,~,~]...
    = AtlasGuidedEM_kubo(XteRkid,atlasRkid,SRkid,K2,Rkidmask,siz,30);

Imap2 = zeros(siz)+4;
Imap2(M3GT == 0) = 0;
Imap2(Imapbla == 1) = 1;  Imap2(Imapbla == 2) = 1;
Imap2(ImapLkid == 1) = 2; Imap2(ImapLkid == 2) = 2;
Imap2(ImapRkid == 1) = 3; Imap2(ImapRkid == 2) = 3;

JI1 = CalcuJI(Imap,M3GT,K2);
disp(JI1);
JI2 = CalcuJI(Imap2,M3GT,K2);
disp(JI2);
%%
imagesc(Imap2(:,:,200)');
axis tight equal off

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
imagesc(Imap2(:,:,200)');
axis tight equal off
caxis([0 4])
colormap(map);
%%
%In = M3E2; InGT = blaGT;
%mu = Sbla.mu; sigma = sqrt(Sbla.Sigma);
%mu = GMMMubla; sigma = sqrt(GMMSigmabla);

%In = M3E2; InGT = LkidGT;
%mu = SLkid.mu; sigma = sqrt(SLkid.Sigma);
%mu = GMMMuLkid; sigma = sqrt(GMMSigmaLkid);

%In = M3E2; InGT = RkidGT;
%mu = SRkid.mu; sigma = sqrt(SRkid.Sigma);
%mu = GMMMuRkid; sigma = sqrt(GMMSigmaRkid);

edge =[0 0:0.01:0.7 0.7];
xlim([0 0.7])

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
edge =[0 0:0.01:0.7 0.7];
In = M3E2; InGT = RkidGT;
hold on
histogram(In(InGT ==1),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'EdgeAlpha',0.4);

%%
blamask2 = zeros(siz); blamaxcomp = zeros(siz);
L1 = bwconncomp(Imap2 == 1);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
blamaxcomp(L1.PixelIdxList{idx}) = 1; 
tmp = bwdist(logical(blamaxcomp)) <  power(bwarea(blamaxcomp(:))/4/pi*3,1/3);
blamask2(tmp) = 1;
blamask2 = logical(and(blamask2,mask3)); 

Lmaxcomp = zeros(siz); Lkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 2);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Lmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Lmaxcomp)) < power(bwarea(Lmaxcomp(:))/4/pi*3,1/3);
Lkidmask2(tmp) = 1;
Lkidmask2 = logical(and(Lkidmask2,mask3));

Rmaxcomp = zeros(siz); Rkidmask2 = zeros(siz);
L1 = bwconncomp(Imap2 == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
Rmaxcomp(L1.PixelIdxList{idx}) = 1;
tmp = bwdist(logical(Rmaxcomp)) <  power(bwarea(Rmaxcomp(:))/4/pi*3,1/3);
Rkidmask2(tmp) = 1;
Rkidmask2 = logical(and(Rkidmask2,mask3));
LRAND = and(Rkidmask2,Lkidmask2);
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
Rkidmask2 = logical(Rkidmask2 - LRAND + tmp2); 

tmp = zeros(siz);
tmp(RRkid > RLkid) = 1;
tmp2 = and(tmp,LRAND);
Lkidmask2 = logical(Lkidmask2 - LRAND + tmp2);

%%
imagesc(M3E1(:,:,200)');
axis tight equal
caxis([0.2 0.4])
%%
%GC bladder
clearvars PPout GraphModel
masktm = blamask;
masktm2 = blamask2;
GT = zeros(siz); GT(M3GT == 1) = 1;
PPorgan = zeros(siz); PPorgan(masktm) = PPbla(:,1) + PPbla(:,2);
PPorgan = imgaussfilt3(PPorgan,5);
PPback = 1.0 - PPorgan; 
PPout(:,1) = PPorgan(masktm2); PPout(:,2) = PPback(masktm2);
PPout = -log(PPout+eps);
PPout = real(PPout);
GraphModel = CreateFullyConnectedGraphWithMask(masktm2);

Kmat = [1,1; 1,1]; K3 = 2;
sumIm = (M3E1 + M3E2)/2;
sumIm = sumIm(masktm2);
%[lambda,h] =ndgrid(100:100:1000 ,0.01:0.1:10);

%lambda = lambda(:);
%h = h(:);
%OutputJI = zeros(size(h,1),1);
%%
%Reaginal term
clearvars PPout
masktm = Lkidmask;
masktm2 = Lkidmask2;
GT = zeros(siz); GT(M3GT == 2) = 1;
PPorgan = zeros(siz); PPorgan(masktm) = ((PPLkid(:,1)+PPLkid(:,2)));
PPorgannew = zeros(siz);
PPorgannew(masktm2) = PPorgan(masktm2);
PPorgan = imgaussfilt3(PPorgannew,5);
PPback = 1.0 - PPorgan;
PPout(:,1) = PPorgan(masktm2); PPout(:,2) = PPback(masktm2);
PPout = -log(PPout+eps);

Kmat = [1,1; 1,1];
K3 = 2;
%sumIm = (M3E1 + M3E2 + M3E3 + M3E4)/4;
sumIm = (M3E2+M3E1)/2;
sumIm = sumIm(masktm2);
%%
imagesc(PPorgan(:,:,200)');

%%
%Reaginal term
clearvars PPout 
masktm = Rkidmask;
masktm2 = Rkidmask2;
GT = zeros(siz); GT(M3GT == 3) = 1;
PPorgan = zeros(siz); PPorgan(masktm) = PPRkid(:,1)+PPRkid(:,2);
PPorgannew = zeros(siz);
PPorgannew(masktm2) = PPorgan(masktm2);
PPorgan = imgaussfilt3(PPorgannew,5);
PPback = 1.0 - PPorgan;
PPout(:,1) = PPorgan(masktm2); PPout(:,2) = PPback(masktm2);
PPout = -log(PPout+eps);

Kmat = [1,1; 1,1];
K3 = 2;
sumIm = (M3E1 + M3E2 )/2;
sumIm = sumIm(masktm2);
%%
se = strel('sphere',5);
I1 = zeros(siz);
I1(Imap2 == 2) = 1;
I1 = imopen(I1,se); I1 = imclose(I1,se);
%%
LK = signed_EUDT(I1); 
LK = LK(mask); 
E1 = sqrt((1 - (LK(GraphModel.Hj)-LK(GraphModel.Hi))./GraphModel.dist)/2);
E2 = sqrt((1 - (RK(GraphModel.Hj)-RK(GraphModel.Hi))./GraphModel.dist)/2);
E1 = real(E1); E2 = real(E2);

%%
imagesc(M3E2(120:220,190:290,200)');
axis tight equal
%caxis([0.04 0.5])
colormap(gray)
%%
[sig,lambda] =ndgrid(0.00005:0.005:0.2,0.2:0.1:0.4);
lambda = lambda(:);
sig = sig(:);
OutputJI = zeros(size(sig,1),1);


%%
shpepri = Iw(masktm2); 
te = (1-(shpepri(GraphModel.Hj)-shpepri(GraphModel.Hi))./GraphModel.dist)./2;
%%
te =real(sqrt(te));
%%
clearvars GraphModel
GraphModel = CreateFullyConnectedGraphWithMask(masktm2);
%%
%�t��
%for n = 1:size(sig,1)
n = 1;
lambda = 0.5;
sig = 0.00505;
    N = size(PPout,1);
    CurLabel = zeros(N,1)+2;
    PreLabel = zeros(N,1);
    Output = zeros(siz);
    PropLabel = ones(N,1);
    
    flag = 0; PreE = 0;
    c = 0.2;
 
    while(flag~=1)
        GraphModel.Vs = PPout((1:N)'+(CurLabel(:)-1)*N)*lambda(n);
        GraphModel.Vt = PPout((1:N)'+(PropLabel(:)-1)*N)*lambda(n);
        
        Z = (sumIm(GraphModel.Hj)-sumIm(GraphModel.Hi)).^2;
        
        %GraphModel.H01(:) = exp(-Z./ (2*sig(n)^2)) ./ GraphModel.dist;
        %GraphModel.H01(:) = c*(exp(-Z./ (2*sig(n)^2)) ./ GraphModel.dist) +(1-c)*te; 
        tem1 = c.*(exp(-Z./ (2*sig(n)^2)) ./ GraphModel.dist); tem2 = (1-c).*te;
        GraphModel.H01(:) = tem1 + tem2;
        GraphModel.H10(:) = GraphModel.H01;
        GraphModel.H00(:) = 0;
        GraphModel.H11(:) = 0;
        [lowerBound, label] = qpboMex([GraphModel.Vs,GraphModel.Vt],...
            [GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]); 
        
        label = logical(label);
        CurLabel(label) = PropLabel(label);
               
          Eunary = sum(GraphModel.Vs); Epairwise1 = sum(tem1); Epairwise2 = sum(tem2);
        E = Eunary + Epairwise1 + Epairwise2;
        disp(E);
        
        if CurLabel == PreLabel
            flag = 1;
            Output(masktm2) = CurLabel;
        end
        
        PreLabel = CurLabel;
        
  
    end
    
    JI3= CalcuJI(Output,GT,1);
    disp(JI3);
    OutputJI(n) = JI3';
    disp(n);
end

%%
save_raw(Output,'C:\\Users\\yourb\\Desktop\\M3GC.raw','*uint8')
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

%%
aaa = zeros(siz); aaa(Lkidmask2) = M3GT(Lkidmask2);
%%
imagesc(Output(:,:,200)');
axis tight equal off

%%
clearvars GraphModel
GraphModel = CreateFullyConnectedGraphWithMask(masktm2);
%�N��
%for n = 1:size(sig,1)
n =1;
lambda = 5;
sig = 0.05;
    N = size(PPout,1);
    CurLabel = zeros(N,1)+2;
    PreLabel = zeros(N,1);
    Output = zeros(siz);
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
        Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H00);
        E = Eunary + Epairwise;
        disp(E);
       
        if CurLabel == PreLabel
            flag = 1;
            Output(masktm2) = CurLabel;
        end
        
        PreLabel = CurLabel;
        
        PreE = E;
    end
    
    JI3= CalcuJI(Output,GT,1);
    disp(JI3);
    OutputJI(n) = JI3';
    disp(n);
%end

%%
Result = zeros(siz) +4;
Result(M3GT == 0 ) = 0;
%%
Result(Output == 1) = 3;
%%
imagesc(Imap2(:,:,70)');
axis tight equal off
caxis([0 4])
colormap(map);
%%
JI4= CalcuJI(Result,M3GT,1);
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];