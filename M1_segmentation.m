
%%
st = 164; en = 439;
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
%initial_value
K=4;
for k = 1:K
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2); tmp3 = Xtr(:,3);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
clearvars tmp1 tmp2 tmp3


%%
%atlas mouse1
sig1 = 5; sig2 = 3; K = 4;
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask1,pM2GT,pM3GT);
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pmask1,siz2);
JI= CalcuJI(Imap,pM1GT,3);
disp(JI);

%%
%voronoi
voronoiIn = zeros(siz2);

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

voronoiOut = zeros(siz2);
[Out,~] = mistVoronoiDistanceTransform(uint8(voronoiIn(pmask1)));
voronoiOut(pmask1) = Out;

%%
GraphModel = CreateFullyConnectedGraphWithMask(pmask1);
%%
I1 = pM1GT ==2;
I1 = logical(I1);

I2 = pM1GT ==3;
I2 = logical(I2);
%%
se = strel('sphere',5);
I1 = zeros(siz2);
I1(pmask1) = L==2;
I1 = imopen(I1,se); I1 = imclose(I1,se);
I1 = logical(I1);

I2 = zeros(siz2);
I2(pmask1) = L ==3;
I2 = imopen(I2,se); I2 = imclose(I2,se);
I2 = logical(I2);
%%
LK = signed_EUDT(I1); RK = signed_EUDT(I2);
LK = LK(pmask1); RK = RK(pmask1);
E1 = sqrt((1 - (LK(GraphModel.Hj)-LK(GraphModel.Hi))./GraphModel.dist)/2);
E2 = sqrt((1 - (RK(GraphModel.Hj)-RK(GraphModel.Hi))./GraphModel.dist)/2);
E1 = real(E1); E2 = real(E2);
%%
test = zeros(siz2);
test(pmask1) = LK;
%%
imagesc(test(:,:,154)');
caxis([-10 10]);
%%
RP = cell(1,K);
for k = 1:K
    RP{1,k} = -log(PP(:,k)+eps);
    RP{1,k}(isnan(RP{1,k})) = 0;
end
RP = cell2mat(RP);
RP = real(RP);

%%
graydiff = zeros(K);
graydiff(1,2) = 1; graydiff(1,3) = 1; graydiff(2,1) = 1; graydiff(3,1) = 1; graydiff(2,3) = 1;graydiff(3,2) = 1;
graydiff(1,4) = 2; graydiff(2,4) = 2; graydiff(3,4) = 2;
graydiff(4,1) = 3; graydiff(4,2) = 3; graydiff(4,3) = 3;

shape = zeros(4);
shape(1,2) = 1; shape(1,3) = 1; shape(2,1) = 1; shape(3,1) = 1; shape(2,3) = 1;shape(3,2) = 1;
shape(1,4) = 0; shape(4,1) = 0;
shape(2,4) = 2; shape(4,2) = 2;
shape(3,4) = 3; shape(4,3) = 3;

%%
[lambda,h,c] =ndgrid(0.01:0.01:0.04,0.1:0.2:1,0.2:0.2:1.0);
lambda = lambda(:);
h = h(:);
c = c(:);
sumJI = zeros(size(h,1),1);

%%
%for n = 1:100
lambda = 0.04;
h = 0.1;
c=1000.0;

N = size(RP,1);
CurLabel = zeros(N,1)+K;
PreLabel = zeros(N,1);
Output = zeros(siz2);
flag = 0;
PreE = 0;
Sigmat =  abs(bsxfun(@minus,GMMMu(:,1),GMMMu(:,1)'))*h + eye(K);
PropLabel = double(Out);
PropLabel(PropLabel == 0) = 1;

while(flag ~=1)
    GraphModel = SetTWeights(GraphModel,RP,CurLabel,PropLabel,lambda,N);
    GraphModel = SetNWeights(GraphModel,pM1E1(pmask1),CurLabel,PropLabel,Sigmat,graydiff,shape,E1,E2,c);
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
%end

JI = zeros(3,1);
for k = 1:K-1
    A = Output == k;   A = A(:);
    B = pM1GT == k;    B = B(:);
    JI(k,1) =  sum(and(A,B)) ./ sum(or(A,B));
end
clearvars A B
disp(JI);

%disp(mean(JI));
sumJI(n) = sum(JI(:));
disp(n);
end

%%
%figure
slice = 191;
subplot(1,4,1);
imagesc(pM1E1(:,:,slice)'); 
axis equal tight off
colormap(gca,'gray');
caxis([0 0.7])

subplot(1,4,2);
imagesc(pM1GT(:,:,slice)'); 
axis equal tight off
colormap(gca,'default')
caxis([0 4]);
%caxis([0 0.05])

subplot(1,4,3);

imagesc(Imap(:,:,slice)'); 
axis equal tight off
colormap(gca,'default')
caxis([0 4]);
%caxis([0 0.05])

subplot(1,4,4);
imagesc(Output(:,:,slice)'); 
axis equal tight off
colormap(gca,'default');
caxis([0 4])

%%
test = pM1GT;
test(test == 4) = 0;

out = Output;
out(out == 4) = 0;

save_raw(pM1E1,'C:\Users\yourb\Desktop\new2\m1.raw','*double');
save_raw(test,'C:\Users\yourb\Desktop\new2\m1gt.raw','*double');
save_raw(out,'C:\Users\yourb\Desktop\new2\outm1.raw','*double');

test1 = zeros(siz2);
test2 = zeros(siz2);
test1(pmask1) = PP(:,2);
test2(pmask1) = PP(:,3);

save_raw(test1,'C:\Users\yourb\Desktop\new2\PP2.raw','*double');
save_raw(test2,'C:\Users\yourb\Desktop\new2\PP3.raw','*double');
%%
seikai = pM1GT(pmask1);
seikai = double(seikai);
%%
Label = CurLabel;
term1 = RP((1:N)'+(Label(:)-1)*N);

Img = pM1E1(pmask1);
Sig = sub2ind(size(Sigmat), Label(GraphModel.Hi), Label(GraphModel.Hj));
Z = Img(GraphModel.Hi)-Img(GraphModel.Hj);
Zval = Z./Sigmat(Sig);
term2 = func(Zval,GraphModel.dist);
term2 = func(-Zval,GraphModel.dist);

shapeval = shape(Sig);
shapeval(shapeval == 2) = E1(shapeval == 2);
shapeval(shapeval == 3) = E2(shapeval == 3);
term3 = shapeval;
%%
term1da = term1(Label == 3);
term2da = term2(Label == 3);
term3da = term3(Label == 3);

disp(sum(term1da(:))*lambda); 
disp(sum(term2da(:))./(1+c));
disp((sum(term3da(:))*c)./(1+c));
term = sum(term1da(:))+sum(term2da(:))+sum(term3da(:));
disp(term);
%%
%%
Label = CurLabel;
E1test = pM1E1(pmask1);
P = GraphModel.Hi(Label(GraphModel.Hi)== 2 & Label(GraphModel.Hj)==4);
Q = GraphModel.Hj(Label(GraphModel.Hi)== 2 & Label(GraphModel.Hj)==4);
%%
sig = 0.0174;
dif = (E1test(P) - E1test(Q))/( sig);
edges = [-2 -2:0.05:3 3];
xlim([-2 3]);
yyaxis left
histogram(dif,edges);
%histogram(dif);

x = edges';
z = exp(-x.^2 /2);
z(x < 0) = 1;

yyaxis right
xlabel('I_p - I_q / \sigma')
plot(x,z,'LineWidth',3);
