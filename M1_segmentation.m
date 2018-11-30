pM1E1 = M1E1(:,:,st:en); pM1E2 = M1E2(:,:,st:en); pM1E3 = M1E3(:,:,st:en); pM1E4 = M1E4(:,:,st:en); pM1GT = M1GT(:,:,st:en); pmask1 = mask1(:,:,st:en);
pM2E1 = M2E1(:,:,st:en); pM2E2 = M2E2(:,:,st:en); pM2E3 = M2E3(:,:,st:en); pM2E4 = M2E4(:,:,st:en); pM2GT = M2GT(:,:,st:en); pmask2 = mask2(:,:,st:en);
pM3E1 = M3E1(:,:,st:en); pM3E2 = M3E2(:,:,st:en); pM3E3 = M3E3(:,:,st:en); pM3E4 = M3E4(:,:,st:en); pM3GT = M3GT(:,:,st:en); pmask3 = mask3(:,:,st:en);
siz2 = size(pM1E1);
clearvars M1E1 M1E2 M1E3 M1E4 M2E1 M2E2 M2E3 M2E4 M3E1 M3E2 M3E3 M3E4 mask1 mask2 mask3 M1GT M2GT M3GT

%%
%train_mouse2_mouse3 test_mouse1
Xtr = [[pM2E2(pmask2); pM3E2(pmask3)] [pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E2(pmask1) pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];
XGTte = [pM1GT(pmask1)];

%%
%train_mouse2_mouse3 test_mouse1
Xtr = [[pM2E2(pmask2); pM3E2(pmask3)] [pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E2(pmask1) pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [maskM2(pmask2); maskM3(pmask3)];
XGTte = [maskM1(pmask1)];

%%
%initial_value
for k = 1:4
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
   % S.mu(k,:) = S.mu(k,:) - 2*sqrt(diag(S.Sigma(:,:,k)))';
    %S.Sigma(:,:,k) = (sqrt(S.Sigma(:,:,k))./2).^2;
end
% S.mu(1,:) = S.mu(1,:) -2*sqrt(diag(S.Sigma(:,:,1)))';
 clearvars tmp1 tmp2 tmp3
%%
%Atlas_guided EM
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask1,pM2GT,pM3GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood] = AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask1,siz2,30);
JI= CalcuJI(Imap,pM1GT,K-1);
disp("EM_MAP result")
disp(JI);
%%
histogram(pM1E2(pM1GT ==1));
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
temp = Imap;
temp(Imap == 4) = 0;

%%
imagesc(temp(:,:,46)');
axis tight equal off
%colormap(gray)
caxis([0 4]);
%%
imagesc(pM1E2(110:300,140:330,206)');
axis tight equal off
colormap(gray)
caxis([0 0.7])
%%
pM1GT(pM1GT ==4) =0;
imagesc(pM1GT(110:300,140:330,206)');
axis tight equal off
colormap(map)
caxis([0 4])

%%
imagesc(voronoiFig(:,:,110)');
axis tight equal off
caxis([0 4])
%%
% for n = 1:27
n = 1;
N = size(RP,1);
CurLabel = zeros(N,1)+K;
PreLabel = zeros(N,1);
Output = zeros(siz2);
flag = 0;
PreE = 0;
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
%subplot(2,1,1)
imagesc(Imapout(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)
%%
subplot(2,1,2)
imagesc(pM1E1(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

