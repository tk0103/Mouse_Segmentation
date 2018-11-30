st = 164; en = 439;
pM1E1 = M1E1(:,:,st:en); pM1E2 = M1E2(:,:,st:en); pM1E3 = M1E3(:,:,st:en); pM1E4 = M1E4(:,:,st:en); pM1GT = M1GT(:,:,st:en); pmask1 = mask1(:,:,st:en);
pM2E1 = M2E1(:,:,st:en); pM2E2 = M2E2(:,:,st:en); pM2E3 = M2E3(:,:,st:en); pM2E4 = M2E4(:,:,st:en); pM2GT = M2GT(:,:,st:en); pmask2 = mask2(:,:,st:en);
pM3E1 = M3E1(:,:,st:en); pM3E2 = M3E2(:,:,st:en); pM3E3 = M3E3(:,:,st:en); pM3E4 = M3E4(:,:,st:en); pM3GT = M3GT(:,:,st:en); pmask3 = mask3(:,:,st:en);
siz2 = size(pM1E1);
%clearvars M1E1 M1E2 M1E3 M1E4 M2E1 M2E2 M2E3 M2E4 M3E1 M3E2 M3E3 M3E4 mask1 mask2 mask3 M1GT M2GT M3GT
%%
%train_mouse1 mouse3 test_mouse2
Xtr = [[pM1E2(pmask1); pM3E2(pmask3)] [pM1E3(pmask1); pM3E3(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E2(pmask2) pM2E3(pmask2) pM2E4(pmask2)];
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
XGTte = [pM2GT(pmask2)];

%%
%train_mouse1 mouse3 test_mouse2
Xtr = [[pM1E2(pmask1); pM3E2(pmask3)] [pM1E3(pmask1); pM3E3(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E2(pmask2) pM2E3(pmask2) pM2E4(pmask2)];
XGTtr = [maskM1(pmask1); maskM3(pmask3)];
XGTte = [maskM2(pmask2)];
%%
imagesc(pM2GT(:,:,115)');
axis tight equal off
%%
%initial_value
K=4;
sig1 = 5; %bladder
sig2 = 3; %kidneys

for k = 1:4
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); 
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
  %  S.Sigma(:,:,k) = (sqrt(S.Sigma(:,:,k))./4).^2;
end
 %S.Sigma(:,:,1) = (sqrt(S.Sigma(:,:,1))./8).^2;
 %S.Sigma(:,:,2) = (sqrt(S.Sigma(:,:,2))./4).^2;
 %S.Sigma(:,:,3) = (sqrt(S.Sigma(:,:,3))./4).^2;
%S.mu(1,:) = S.mu(1,:) +1*sqrt(diag(S.Sigma(:,:,1)))';
%S.mu(2,:) = S.mu(2,:) +1*sqrt(diag(S.Sigma(:,:,2)))';
%S.mu(3,:) = S.mu(3,:) +1*sqrt(diag(S.Sigma(:,:,3)))';

clearvars tmp1 tmp2 tmp3
%%
%Atlas_guided EM2
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask2,pM1GT,pM3GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood] = ...
    AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask2,siz2,30);
JI= CalcuJI(Imap,pM2GT,K-1);
disp("EM_MAP result")
disp(JI);
%%
temp = zeros(siz2);
temp(pmask2) = atlas(:,1);
%%
temp = Imap;
temp(Imap == 4) = 0;
%%
imagesc(Imap2(:,:,271)');
axis tight equal off
%caxis([0 0.00000000001])
%colormap(gray)
%caxis([0 0.7])

%%
%Reaginal term
tempPP1 = zeros(siz2); tempPP1(mask) = PP(:,1); 
tempPP1 = imgaussfilt3(tempPP1,5);
tempPP2 = imgaussfilt3(PPtemp1,10);
tempPP3 = imgaussfilt3(PPtemp2,10);
tempPP4 = 1 - tempPP1 - tempPP2 - tempPP3;

PPout(:,1) = tempPP1(mask); PPout(:,2) = tempPP2(mask); 
PPout(:,3) = tempPP3(mask); PPout(:,4) = tempPP4(mask);
%clearvars tempPP1 tempPP2 tempPP3 tempPP4
%%
temp = zeros(siz2);
temp(pmask2) = PP(:,1);
%%
slice1 = 180;
slice2 = 66;
imagesc(tempPP4(:,:,slice1)');
axis tight equal off
colormap(map)
%caxis([0 4]);



%%
RP = cell(1,K);
for k = 1:K
    RP{1,k} = -log(PPout(:,k)+eps);
   % RP{1,k}(isnan(RP{1,k})) = 0;
end
RP = cell2mat(RP);
%RP = real(RP);
%%
temp = zeros(siz2);
temp(pmask2) = RP(:,4);
%%
imagesc(temp(:,:,216)');
axis tight equal off

%%
%GraphCut
GraphModel = CreateFullyConnectedGraphWithMask(pmask2);

%Shape term
[E1,E2] = Create_shape(L,GraphModel,pmask2,siz2);

%voronoi
voronoiIn = zeros(siz2);
for k = 1:K-1
    LL = bwconncomp(Imap2 == k);
    numPixels = cellfun(@numel,LL.PixelIdxList);
    [~,idx] = max(numPixels);
    voronoiIn(LL.PixelIdxList{idx}) = k;
end
voronoiFig = zeros(siz2);
[voronoiOut,~] = mistVoronoiDistanceTransform(uint8(voronoiIn(pmask2)));
voronoiFig(pmask2) = voronoiOut;

%%
%for n =1:12
n = 1;
lambda = 0.1;
h =1;
c = 0;
N = size(RP,1);
CurLabel = zeros(N,1)+K;
PreLabel = zeros(N,1);
Output = zeros(siz2);
flag = 0;
PreE = 0;
Sigmat =  abs(bsxfun(@minus,GMMMu(1:K,1),GMMMu(1:K,1)'))*h(n) + eye(K);
PropLabel = double(voronoiOut);
PropLabel(PropLabel == 0) = 1;
Kmat = ones(4);
while(flag ~=1)
    %for aa = 1:4
    %    PropLabel = zeros(N,1) + aa;
    GraphModel = SetTWeights(GraphModel,RP,CurLabel,PropLabel,lambda(n),N);
  %  GraphModel = SetNWeights(GraphModel,pM2E1(pmask2),CurLabel,PropLabel,Sigmat,graydiff,shape,E1,E2,c(n));
    GraphModel = SetNWeights(GraphModel,pM2E1(mask),CurLabel,PropLabel,Sigmat,Kmat);
     
   [lowerBound, labels] = qpboMex([GraphModel.Vs,GraphModel.Vt],[GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
    labels = logical(labels);
    CurLabel(labels) = PropLabel(labels);
    
 %   GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
   % Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H00);
   % E = Eunary + Epairwise;
   % end 
    
    disp(lowerBound);
    if CurLabel == PreLabel
        flag = 1;
        Output(mask) = CurLabel;
    end
    
    PreLabel = CurLabel;
    PreE = lowerBound;
end

JI= CalcuJI(Output,pM2GT,K-1);
disp("GraphCut_JI")
disp(JI);

OutputJI(n,:) = JI';
disp(n);
%end

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
%%
Imapout = Output;
Imapout(Output==4) = 0;
pM2GTout = pM2GT;
pM2GTout(pM2GT==4) = 0;
%%
save_raw(Output,'C:\\Users\\yourb\\Desktop\\M2GCwork.raw','*uint8')
%%
slice1 = 226;
slice2 = 66;

subplot(2,2,1)
imagesc(Imapout(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,2)
imagesc(voronoiFig(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,3)
imagesc(Imapout(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,4)
imagesc(voronoiFig(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

%%
subplot(2,1,1)
imagesc(pM2E1(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

subplot(2,1,2)
imagesc(pM2E1(:,:,slice2)');
axis tight equal off
caxis([0 0.7])
colormap(gray)


%%
imagesc(Imapout(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)