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
XGTte = pM3GT(pmask3);
%%
sig1 = 5; sig2 = 3; K = 4;
%initial_value
for k = 1:K
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
S.Sigma(:,:,1) = sqrt(S.Sigma(:,:,1));
clearvars tmp1 tmp2 tmp3

%%
%initial_value test
for k = 1:K
    tmp1 = Xte(:,1); tmp2 = Xte(:,2); tmp3 = Xte(:,3);
    STe.mu(k,1) = mean(tmp1(XGTte == k));
    STe.mu(k,2) = mean(tmp2(XGTte == k));
    STe.mu(k,3) = mean(tmp3(XGTte == k));
    STe.Sigma(:,:,k) = cov(([tmp1(XGTte == k),tmp2(XGTte == k),tmp3(XGTte == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
%Atlas_guided EM
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask3,pM1GT,pM2GT);
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pmask3,siz2);
JI= CalcuJI(Imap,pM3GT,K-1);
disp("EM_MAP result")
disp(JI);
%%
temp = zeros(siz2);
temp(pmask3) = PP(:,1);
%%
imagesc(temp(:,:,slice2)');
axis tight equal off
colormap(gray)%caxis([50 250])
%%
temp1 = pM1GT==1;
temp2 = pM1GT==2;
temp3 = pM1GT==3;
%%
Imap(Imap==4) =0;
%%
imagesc(Imap(:,:,slice1)');
colormap(map)
caxis([0 4])
axis tight equal off
 rectangle('Position',[110,140,190,190],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
%%
imagesc(Imap(110:300,140:330,slice1)');
axis tight equal off
%colormap(gray)
%caxis([0 0.7])
colormap(map)
caxis([0 4])

%%
temp = pM1GT; temp(pM1GT==4) =0;
%%
imagesc(pwM1E1(:,:,slice2)');
colormap(gray)
caxis([0 0.7])
axis tight equal off
 rectangle('Position',[315,250,70,70],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
%%
hold on
contour(temp1(:,:,slice2)',[1 1],'red');
axis tight equal off

%%
hold on
contour(temp2(:,:,slice1)',[1 1],'red');
axis tight equal off

contour(temp3(:,:,slice1)',[1 1],'red');
axis tight equal off

%%
imagesc(pM3E2(100:290,140:330,slice1)');
colormap(gray)
caxis([0 0.7])
axis tight equal off

hold on
contour(temp2(100:290,140:330,slice1)',[1 1],'red');
axis tight equal off

contour(temp3(100:290,140:330,slice1)',[1 1],'red');
axis tight equal off


%%
%Reaginal term
tempPP1 = zeros(siz2); tempPP2 = zeros(siz2); tempPP3 = zeros(siz2);
tempPP1(pmask3) = PP(:,1); tempPP2(pmask3) = PP(:,2); tempPP3(pmask3) = PP(:,3);
tempPP1 = imgaussfilt3(tempPP1,5);
tempPP2 = imgaussfilt3(tempPP2,5);
tempPP3 = imgaussfilt3(tempPP3,5);
tempPP4 = 1 - tempPP1 - tempPP2 - tempPP3;
PPout(:,1) = tempPP1(pmask3); PPout(:,2) = tempPP2(pmask3); 
PPout(:,3) = tempPP3(pmask3); PPout(:,4) = tempPP4(pmask3);
%%
temp1 = pM3E2(245:375,225:355,66);
temp2 = pM3E2(246:376,226:356,66);
temp3 = pM3GT(245:375,225:355,66);
temp4 = zeros(size(temp1)); temp4(temp3 ==1) = 1;
tempmask = pmask3(245:375,225:355,66);
Z = temp1 - temp2;
Zdash = Z.^2;
%%
temp = exp(-(Z.^2) /2);
%temp(Z<0) = 1;
%%
tempout = zeros(size(temp));
tempout(tempmask) =temp(tempmask);

%%
temp = zeros(siz2);
temp(pmask3) = PP(:,1);
temp =  temp(245:375,225:355,66);
%%

imagesc(temp');
axis tight equal off
caxis([0.995 1.0])
colormap(gray)

hold on
contour(temp4',[1 1],'red');
axis tight equal off

%%
RP = cell(1,K);
for k = 1:K
    RP{1,k} = -log(PPout(:,k)+eps);
    RP{1,k}(isnan(RP{1,k})) = 0;
end
RP = cell2mat(RP);
RP = real(RP);


%GraphCut
GraphModel = CreateFullyConnectedGraphWithMask(pmask3);

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
for n =2:125
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
    GraphModel = SetNWeights(GraphModel,pM3E2(pmask3),CurLabel,PropLabel,Sigmat,graydiff,shape,E1,E2,c(n));
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

OutputJI(n,:) = JI';
disp(n);
end
%%
imagesc(Output(:,:,220)');
axis tight equal
%%
save_raw(Imap,'C:\Users\yourb\Desktop\Imap.raw','*uint8');
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
%%
outGC = Output;
outGC(Output==4) = 0;
pM3GTout = pM3GT;
pM3GTout(pM3GT==4) = 0;
%%
slice1 = 206;
slice2 = 66;

subplot(2,2,1)
imagesc(outGC(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,2)
imagesc(pM3GTout(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,3)
imagesc(outGC(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,4)
imagesc(pM3GTout(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

%%
imagesc(outGC(120:250,210:340,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)
%%
imagesc(pM1E2(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
 rectangle('Position',[120,210,130,130],'FaceColor','none','EdgeColor','r',...
    'LineWidth',2)

%%
%subplot(2,1,1)
val = slice1;
imagesc(Imapout(:,:,val)');
axis tight equal off
caxis([0 4])
colormap(map)