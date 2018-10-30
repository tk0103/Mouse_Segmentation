%%
st = 164; en = 439;
pwM1E1 = wM1E1(:,:,st:en); pwM2E1 = wM2E1(:,:,st:en); pwM3E1 = wM3E1(:,:,st:en);
pwM1GT = wM1GT(:,:,st:en);  pwM2GT = wM2GT(:,:,st:en);  pwM3GT = wM3GT(:,:,st:en);
pwmask1 = wmask1(:,:,st:en);  pwmask2 = wmask2(:,:,st:en); pwmask3 = wmask3(:,:,st:en);
siz = size(pwM2E1);
clearvars wM1E1 wM2E1 wM3E1 wM1GT wM2GT wM3GT wmask1 wmask2 wmask3
%%
%train_mouse2_mouse3 test_mouse1
Xtr = [pwM2E1(pwmask2); pwM3E1(pwmask3) ];
Xte = pwM1E1(pwmask1);
XGTtr = [pwM2GT(pwmask2); pwM3GT(pwmask3)];

%%
%train_mouse1 mouse3 test_mouse2
Xtr = [pwM1E1(pwmask1); pwM3E1(pwmask3) ];
Xte = pwM2E1(pwmask2);
XGTtr = [pwM1GT(pwmask1); pwM3GT(pwmask3)];

%%
%train_mouse1 mouse2 test_mouse3
Xtr = [pwM1E1(pwmask1); pwM2E1(pwmask2) ];
Xte = pwM3E1(pwmask3);
XGTtr = [pwM1GT(pwmask1); pwM2GT(pwmask2)];

%%
%initial_value
for k = 1:K
    tmp1 = Xtr(:,1);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.Sigma(:,:,k) = cov(tmp1(XGTtr == k));
end

%%
%Atlas_guided EM,Mouse1
sig1 = 5; sig2 = 3; K = 4;
atlas  = atlasfunc2(sig1,sig2,K,siz,pwmask1,pwM2GT,pwM3GT);
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pwmask1,siz);
JI= CalcuJI(Imap,pwM1GT,K-1);
disp("EM_MAP result")
disp(JI);

%%
%Atlas_guided EM,Mouse2
sig1 = 5; sig2 = 3; K = 4;
atlas  = atlasfunc2(sig1,sig2,K,siz,pwmask2,pwM1GT,pwM3GT);
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pwmask2,siz);
JI= CalcuJI(Imap,pwM2GT,K-1);
disp("EM_MAP result")
disp(JI);
%%
%Atlas_guided EM,Mouse3
sig1 = 5; sig2 = 3; K = 4;
atlas  = atlasfunc2(sig1,sig2,K,siz,pwmask3,pwM1GT,pwM2GT);
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pwmask3,siz);
JI= CalcuJI(Imap,pwM3GT,K-1);
disp("EM_MAP result")
disp(JI);

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
%%
Imapout = Imap;
Imapout(Imap==4) = 0;
%%
pM1GTout = pM3GT;
pM1GTout(pM1GT==4) = 0;
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
subplot(1,2,1)
imagesc(pwM3E1(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

subplot(1,2,2)
imagesc(pwM3E1(:,:,slice2)');
axis tight equal off
caxis([0 0.7])
colormap(gray)


%%
subplot(1,2,1)
imagesc(pM1E2(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

subplot(1,2,2)
imagesc(pM1E2(:,:,slice2)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

%%
imagesc(Imapout(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)
 rectangle('Position',[270,210,110,110],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
%%
imagesc(Imapout(270:380,210:320,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)




%%
%EM initialvalue test
K = 4;
for k = 1:K
    tmp1 = Xte(:,1); tmp2 = Xte(:,2); tmp3 = Xte(:,3); 
    Ste.Mu(k,1) = mean(tmp1(XGTte == k));
    Ste.Mu(k,2) = mean(tmp2(XGTte == k));
    Ste.Mu(k,3) = mean(tmp3(XGTte == k));
    Ste.Sigma(:,:,k) = cov([tmp1(XGTte == k),tmp2(XGTte == k),tmp3(XGTte == k)]);
end

%%
voronoiIn = zeros(siz);
LL = bwconncomp(Imap==1);
numPixels = cellfun(@numel,LL.PixelIdxList);
[~,idx] = max(numPixels);
voronoiIn(LL.PixelIdxList{idx}) = 1;

LL = bwconncomp(Imap==2);
numPixels = cellfun(@numel,LL.PixelIdxList);
[~,idx] = max(numPixels);
voronoiIn(LL.PixelIdxList{idx}) = 2;

LL = bwconncomp(Imap==3);
numPixels = cellfun(@numel,LL.PixelIdxList);
[~,idx] = max(numPixels);
voronoiIn(LL.PixelIdxList{idx}) = 3;

voronoiOut = zeros(siz);
[Out,~] = mistVoronoiDistanceTransform(uint8(voronoiIn(pwmask2)));
voronoiOut(pwmask2) = Out;

%%
[lambda,h,c] =ndgrid(0.05:0.25:0.8,0.1:0.2:1,0.2:0.2:2.0);
lambda = lambda(:);
h = h(:);
c = c(:);
sumJI = zeros(size(h,1),1);
%%
GraphModel = CreateFullyConnectedGraphWithMask(pwmask2);

%%
I1 = zeros(siz);
I1(pwmask3) = L==2;
se = strel('sphere',5);
I1 = imopen(I1,se);
I1 = imclose(I1,se);
I1 = logical(I1);

I2 = zeros(siz);
I2(pwmask3) = L==3;
se = strel('sphere',5);
I2 = imopen(I2,se);
I2 = imclose(I2,se);
I2 = logical(I2);

M3LK = signed_EUDT(I1);
M3RK = signed_EUDT(I2);

M3LKIn = M3LK(pwmask3);
M3RKIn = M3RK(pwmask3);
E1 = sqrt((1 - (M3LKIn(GraphModel.Hj)-M3LKIn(GraphModel.Hi))./GraphModel.dist)/2);
E2 = sqrt((1 - (M3RKIn(GraphModel.Hj)-M3RKIn(GraphModel.Hi))./GraphModel.dist)/2);
E1 = real(E1);
E2 = real(E2);
%%
imagesc(M3LK(:,:,245)');
caxis([0 20]);
%%
RP = cell(1,K);
for k = 1:K
    RP{1,k} = -log(PP(:,k)+eps);
    RP{1,k}(isnan(RP{1,k}))=0;
end
RP = cell2mat(RP);
RP = real(RP);

%%
graydiff = zeros(K);
graydiff(1,2) = 1; graydiff(1,3) = 1; graydiff(2,1) = 1; graydiff(3,1) = 1; graydiff(2,3) = 1;graydiff(3,2) = 1;
graydiff(1,4) = 2; graydiff(2,4) = 2; graydiff(3,4) = 2;
graydiff(4,1) = 3; graydiff(4,2) = 3; graydiff(4,3) = 3;
Sigmat =  abs(bsxfun(@minus,S.mu(:,1),S.mu(:,1)'))*h(n) +eye(K);

shape = zeros(4);
shape(1,2) = 1; shape(1,3) = 1; shape(2,1) = 1; shape(3,1) = 1; shape(2,3) = 1;shape(3,2) = 1;
shape(1,4) = 0; shape(4,1) = 0;
shape(2,4) = 2; shape(4,2) = 2;
shape(3,4) = 3; shape(4,3) = 3;

%%
n=22;
c=1;
%for n = 1:3
    N = size(RP,1);
    CurLabel = zeros(N,1)+K;
    PreLabel = zeros(N,1);
    Output = zeros(siz);
    flag = 0;
    PreE = 0;
    Sigmat =  abs(bsxfun(@minus,GMMMu(:,1),GMMMu(:,1)'))*h(n) + eye(K);
    PropLabel = double(Out);
    PropLabel(PropLabel == 0) = 1;

    while(flag ~=1)
        GraphModel = SetTWeights(GraphModel,RP,CurLabel,PropLabel,lambda(n),N);
        GraphModel = SetNWeights(GraphModel,pwM2E1(pwmask2),CurLabel,PropLabel,Sigmat,graydiff,shape,E1,E2,c);
      
        [lowerBound, labels] = qpboMex([GraphModel.Vs,GraphModel.Vt],[GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
        labels = logical(labels);
        CurLabel(labels) = PropLabel(labels);
        
        GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
        Eunary = sum(GraphModel.Vs);
        Epairwise = sum(GraphModel.H00);
        E = Eunary + Epairwise;
        
        disp(E);
        if CurLabel == PreLabel
            flag = 1;
            Output(pwmask2) = CurLabel;
        end
        
        PreLabel = CurLabel;
        PreE = E;
    end
%%
JI = zeros(1,1);     
for k = 1:3
    A = Output == k;   A = A(:);
    B = pwM3GT == k;    B = B(:);
    JI(k,1) =  sum(and(A,B)) ./ sum(or(A,B));
end
clearvars A B
disp(JI);
%disp(mean(JI));
sumJI(n) = JI;
disp(n);
%end

%%
%1ŽŸŒ³
ori = pwM3E1;
M1B = ori(pwM3GT==1);
M1K = ori(or(pwM3GT==2,pwM3GT==3));
M1E = ori(pwM3GT==4);
M1min = min(ori(:));
M1max = max(ori(:));

h = (M1max - M1min)/ (sqrt(size(ori(pwmask3),1))/5);
edges = M1min:h:M1max;
[NM1B,~] = histcn(M1B,edges);
[NM1K,~] = histcn(M1K,edges);
[NM1E,~] = histcn(M1E,edges);
NM1B = NM1B ./ sum(NM1B) + eps;
NM1K = NM1K ./ sum(NM1K) + eps;
NM1E = NM1E ./ sum(NM1E) + eps;

P = NM1B;
Q = NM1K;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1B;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1K;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);

%%
slice = 206;
y =290;
rangex = [160,230];
rangey = [y-35,y+35];
test1 = pM1E2(:,y,slice);
test2 = pwM1E1(:,y,slice);


subplot(3,4,[1,6])
imagesc(pM1E2(:,:,slice)');
caxis([0,0.7])
colormap gray
axis equal tight 
rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);


subplot(3,4,[3,8])
imagesc(pwM1E1(:,:,slice)');
colormap gray
caxis([0,0.7])
axis equal tight 
 rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);


M1 = movmean(test1,3);
subplot(3,4,[9,10])
h= plot([test1]);
h(1).Color = 'g';
xlim(rangex);
ylim([0,0.4]);
legend({'Gray value'});


M2 = movmean(test2,3);
subplot(3,4,[11,12])
h= plot([test2]);
h(1).Color = 'g';
xlim(rangex);
ylim([0,0.4]);
legend({'Gray value'});

%%
slice = 206;
y =190;
rangex = [175,245];
rangey = [y-35,y+35];
test1 = pM1E2(:,y,slice);
test2 = pwM1E1(:,y,slice);


subplot(3,4,[1,6])
imagesc(pM1E2(:,:,slice)');
caxis([0,0.6])
colormap gray
axis equal tight
rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);


subplot(3,4,[3,8])
imagesc(pwM1E1(:,:,slice)');
colormap gray
caxis([0,0.6])
axis equal tight 
 rectangle('Position',[1,y-0.5,436,1],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
xlim(rangex);
ylim(rangey);


M1 = movmean(test1,3);
subplot(3,4,[9,10])
h= plot([test1]);
h(1).Color = 'g';
xlim(rangex);
ylim([0,0.6]);
legend({'Gray value'});


M2 = movmean(test2,3);
subplot(3,4,[11,12])
h= plot([test2]);
h(1).Color = 'g';
xlim(rangex);
ylim([0,0.6]);
legend({'Gray value'});
%%
x = [0.602,0.637];
y = [3.516,13.33];
scatter(x,y,'.')
legend("a","b")
