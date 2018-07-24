st = 164; en = 439;
pM1E1 = M1E1(:,:,st:en); pM1E2 = M1E2(:,:,st:en); pM1E3 = M1E3(:,:,st:en); pM1E4 = M1E4(:,:,st:en); pM1GT = M1GT(:,:,st:en); pmask1 = mask1(:,:,st:en);
pM2E1 = M2E1(:,:,st:en); pM2E2 = M2E2(:,:,st:en); pM2E3 = M2E3(:,:,st:en); pM2E4 = M2E4(:,:,st:en); pM2GT = M2GT(:,:,st:en); pmask2 = mask2(:,:,st:en);
pM3E1 = M3E1(:,:,st:en); pM3E2 = M3E2(:,:,st:en); pM3E3 = M3E3(:,:,st:en); pM3E4 = M3E4(:,:,st:en); pM3GT = M3GT(:,:,st:en); pmask3 = mask3(:,:,st:en);
siz2 = size(pM1E1);
clearvars M1E1 M1E2 M1E3 M1E4 M2E1 M2E2 M2E3 M2E4 M3E1 M3E2 M3E3 M3E4 mask1 mask2 mask3 M1GT M2GT M3GT

%%
%train_mouse1 mouse2 test_mouse3
Xtr = [[pM1E2(pmask1); pM2E2(pmask2)] [pM1E3(pmask1); pM2E3(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)] ];
Xte = [pM3E2(pmask3) pM3E3(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];

%%
%initial_value
for k = 1:K
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
clearvars tmp1 tmp2 tmp3

%Atlas_guided EM
sig1 = 5; sig2 = 3; K = 4;
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask3,pM1GT,pM2GT);
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pmask3,siz2);
JI= CalcuJI(Imap,pM3GT,K-1);
disp("EM_MAP result")
disp(JI);
%%
%GraphCut
GraphModel = CreateFullyConnectedGraphWithMask(pmask3);

%Reaginal term
RP = cell(1,K);
for k = 1:K
    RP{1,k} = -log(PP(:,k)+eps);
    RP{1,k}(isnan(RP{1,k})) = 0;
end
RP = cell2mat(RP);
RP = real(RP);

%Shape term
[E1,E2] = Create_shape(L,GraphModel,pmask3,siz2);

%voronoi
voronoiIn = zeros(siz2);
for k = 1:K-1
    LL = bwconncomp(Imap==k);
    numPixels = cellfun(@numel,LL.PixelIdxList);
    [~,idx] = max(numPixels);
    voronoiIn(LL.PixelIdxList{idx}) = k;
end
voronoiFig = zeros(siz2);
[voronoiOut,~] = mistVoronoiDistanceTransform(uint8(voronoiIn(pmask3)));
voronoiFig(pmask3) = voronoiOut;

%%
N = size(RP,1);
CurLabel = zeros(N,1)+K;
PreLabel = zeros(N,1);
Output = zeros(siz2);
flag = 0;
PreE = 0;
Sigmat =  abs(bsxfun(@minus,GMMMu(:,1),GMMMu(:,1)'))*h + eye(K);
PropLabel = double(voronoiOut);
PropLabel(PropLabel == 0) = 1;

while(flag ~=1)
    GraphModel = SetTWeights(GraphModel,RP,CurLabel,PropLabel,lambda,N);
    GraphModel = SetNWeights(GraphModel,pM3E1(pmask3),CurLabel,PropLabel,Sigmat,graydiff,shape,E1,E2,c);
    [lowerBound, labels] = qpboMex([GraphModel.Vs,GraphModel.Vt],[GraphModel.Hi,GraphModel.Hj,GraphModel.H00,GraphModel.H01,GraphModel.H10,GraphModel.H11]);
    labels = logical(labels);
    CurLabel(labels) = PropLabel(labels);
    
    GraphModel.Vs(~isfinite(GraphModel.Vs)) = 0;
    Eunary = sum(GraphModel.Vs); Epairwise = sum(GraphModel.H00);
    E = Eunary + Epairwise;
    
    disp(E);
    if CurLabel == PreLabel
        flag = 1;
        Output(pmask3) = CurLabel;
    end
    
    PreLabel = CurLabel;
    PreE = E;
end

JI= CalcuJI(Output,pM3GT,K-1);
disp("GraphCut_JI")
disp(JI);

%sumJI(n) = sum(JI(:));
%disp(n);
%end
%%
% 0.1 0.04
temp = pM3E2-pM3E1;
temp = temp(pmask3);
temp = temp>0.03;

test = zeros(siz2);
test(pmask3) = temp;
%%
imagesc(temtem(:,:,215));
caxis([0 2])
%%
imagesc(GT1(:,:,215));
%%
GT1 = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\M3GTnew.raw','*uint8');
GT1 = reshape(GT1,siz);
GT1 = GT1(:,:,st:en);
GT1 = double(GT1);
%%
temtem = test+GT1;
%%
save_raw(out,'C:\\Users\\yourb\\Desktop\\NZ_unet\\M3GTnewnew.raw')


%%
out = zeros(siz);
out(:,:,st:en) = temtem;
