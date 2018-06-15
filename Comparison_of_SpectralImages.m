pM1E1 = M1E1(:,:,st:en); pM1E2 = M1E2(:,:,st:en); pM1E3 = M1E3(:,:,st:en); pM1E4 = M1E4(:,:,st:en); pM1GT = M1GT(:,:,st:en); pmask1 = mask1(:,:,st:en);
pM2E1 = M2E1(:,:,st:en); pM2E2 = M2E2(:,:,st:en); pM2E3 = M2E3(:,:,st:en); pM2E4 = M2E4(:,:,st:en); pM2GT = M2GT(:,:,st:en); pmask2 = mask2(:,:,st:en);
pM3E1 = M3E1(:,:,st:en); pM3E2 = M3E2(:,:,st:en); pM3E3 = M3E3(:,:,st:en); pM3E4 = M3E4(:,:,st:en); pM3GT = M3GT(:,:,st:en); pmask3 = mask3(:,:,st:en);
siz2 = size(pM1E1);
clearvars M1E1 M1E2 M1E3 M1E4 M2E1 M2E2 M2E3 M2E4 M3E1 M3E2 M3E3 M3E4 mask1 mask2 mask3 M1GT M2GT M3GT

%%
%train_mouse2_mouse3 test_mouse1
Xtr = [pM2E4(pmask2); pM3E4(pmask3)];
Xte = pM1E4(pmask1);
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];

%%
Xtr = [[pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];
%%
Xtr = [[pM2E1(pmask2); pM3E1(pmask3)] [pM2E2(pmask2); pM3E2(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E1(pmask1) pM1E2(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];
%%
Xtr = [[pM2E1(pmask2); pM3E1(pmask3)] [pM2E2(pmask2); pM3E2(pmask3)] [pM2E3(pmask2); pM3E3(pmask3)] [pM2E4(pmask2); pM3E4(pmask3)] ];
Xte = [pM1E1(pmask1) pM1E2(pmask1) pM1E3(pmask1) pM1E4(pmask1)];
XGTtr = [pM2GT(pmask2); pM3GT(pmask3)];

%%
%train_mouse1 mouse3 test_mouse2
Xtr = [pM1E4(pmask1); pM3E4(pmask3)];
Xte = pM2E4(pmask2);
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
%%
Xtr = [[pM1E3(pmask1); pM3E3(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E3(pmask2) pM2E4(pmask2)];
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
%%
Xtr = [[pM1E1(pmask1); pM3E1(pmask3)] [pM1E2(pmask1); pM3E2(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E1(pmask2) pM2E2(pmask2) pM2E4(pmask2)];
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
%%
Xtr = [[pM1E1(pmask1); pM3E1(pmask3)] [pM1E2(pmask1); pM3E2(pmask3)] [pM1E3(pmask1); pM3E3(pmask3)] [pM1E4(pmask1); pM3E4(pmask3)] ];
Xte = [pM2E1(pmask2) pM2E2(pmask2) pM2E3(pmask2) pM2E4(pmask2)];
XGTtr = [pM1GT(pmask1); pM3GT(pmask3)];
%%
%train_mouse1 mouse2 test_mouse3
Xtr = [pM1E4(pmask1); pM2E4(pmask2)];
Xte = pM3E4(pmask3);
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];
%%
Xtr = [[pM1E3(pmask1); pM2E3(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)]];
Xte = [pM3E3(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];
%%
Xtr = [[pM1E1(pmask1); pM2E1(pmask2)] [pM1E2(pmask1); pM2E2(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)] ];
Xte = [pM3E1(pmask3) pM3E2(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];
XGTte = [pM3GT(pmask3)];
%%
Xtr = [[pM1E1(pmask1); pM2E1(pmask2)] [pM1E2(pmask1); pM2E2(pmask2)] [pM1E3(pmask1); pM2E3(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)] ];
Xte = [pM3E1(pmask3) pM3E2(pmask3) pM3E3(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];


%%
%initial_value
for k = 1:K
    tmp1 = Xtr(:,1);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.Sigma(:,:,k) = cov((tmp1(XGTtr == k)));
    S.ComponentProportion(k,1) = numel(tmp1(XGTtr == k));
end
S.ComponentProportion =  S.ComponentProportion ./(sum( S.ComponentProportion));
%%
%initial_value
for k = 1:K
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k)]));
    S.ComponentProportion(k,1) = numel(tmp1(XGTtr == k));
end
S.ComponentProportion =  S.ComponentProportion ./(sum( S.ComponentProportion));
%%
for k = 1:K
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
    S.ComponentProportion(k,1) = numel(tmp1(XGTtr == k));
end
S.ComponentProportion =  S.ComponentProportion ./(sum( S.ComponentProportion));
%%
for k = 1:K
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); tmp4 = Xtr(:,4);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.mu(k,4) = mean(tmp4(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k),tmp4(XGTtr == k)]));
    S.ComponentProportion(k,1) = numel(tmp1(XGTtr == k));
end
S.ComponentProportion =  S.ComponentProportion ./(sum( S.ComponentProportion));

%%
%atlas mouse1
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask1,pM2GT,pM3GT);
%%
%atlas mouse2
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask2,pM1GT,pM3GT);

%%
%atlas mouse3
atlas  = atlasfunc2(sig1,sig2,K,siz2,pmask3,pM1GT,pM2GT);

%%
%EM_MAP
GMModel = fitgmdist(Xte,K,'start',S,'ProbabilityTolerance',1e-6,'Options',statset('Display','iter','MaxIter',15));
Mu = GMModel.mu;
Sig = GMModel.Sigma;
proportion = GMModel.ComponentProportion;
y1 = proportion(1) .* mvnpdf(Xte,Mu(1,:),Sig(:,:,1));
y2 = proportion(2) .* mvnpdf(Xte,Mu(2,:),Sig(:,:,2));
y3 = proportion(3) .* mvnpdf(Xte,Mu(3,:),Sig(:,:,3));
y4 = proportion(4) .* mvnpdf(Xte,Mu(4,:),Sig(:,:,4));
F = [y1,y2,y3,y4];
p_l = atlas;
PP = bsxfun(@times,F,p_l);
p_x = sum(PP,2);
r = p_x>0;
PP(r,:) = bsxfun(@rdivide,PP(r,:),p_x(r,:));
PP = bsxfun(@rdivide,PP,p_x);
[~,L] = max(PP,[],2);
I = zeros(siz2);
I(pmask3) = L;

%%
%Atlas_guided EM Mouse1
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pmask1,siz2);
JI= CalcuJI(Imap,pM1GT,K-1);
disp("AtlasGuided EM_MAP result"); disp(JI);

%%
%Atlas_guided EM Mouse2
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pmask2,siz2);
JI= CalcuJI(Imap,pM2GT,K-1);
disp("AtlasGuided EM_MAP result"); disp(JI);
%%
%Atlas_guided EM Mouse3
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,pmask3,siz2);
JI= CalcuJI(Imap,pM3GT,K-1);
disp("AtlasGuided EM_MAP result"); disp(JI);