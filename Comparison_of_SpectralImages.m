%Train M2 M3 M4 Test M1
%Number of Energy 1
Xtr = [M2E1(mask2); M3E1(mask3); M4E1(mask4)];
Xte = M1E1(mask1);
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M2E2(mask2); M3E2(mask3); M4E2(mask4)];
Xte = M1E2(mask1);
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M2E3(mask2); M3E3(mask3); M4E3(mask4)];
Xte = M1E3(mask1);
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M2E4(mask2); M3E4(mask3); M4E4(mask4)];
Xte = M1E4(mask1);
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E2(mask2); M3E2(mask3); M4E2(mask4)]];
Xte = [M1E1(mask1) M1E2(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E3(mask2); M3E3(mask3); M4E3(mask4)]];
Xte = [M1E1(mask1) M1E3(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E1(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M2E2(mask2); M3E2(mask3); M4E2(mask4)] [M2E3(mask2); M3E3(mask3); M4E3(mask4)]];
Xte = [M1E2(mask1) M1E3(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M2E2(mask2); M3E2(mask3); M4E2(mask4)] [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E2(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M2E3(mask2); M3E3(mask3); M4E3(mask4)] [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E3(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E2(mask2); M3E2(mask3);M4E2(mask4)] ...
      [M2E3(mask2); M3E3(mask3); M4E3(mask4)]];
Xte = [M1E1(mask1) M1E2(mask1) M1E3(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E2(mask2); M3E2(mask3);M4E2(mask4)] ...
      [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E1(mask1) M1E2(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M2E2(mask2); M3E2(mask3); M4E2(mask4)] [M2E3(mask2); M3E3(mask3);M4E3(mask4)] ...
      [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E2(mask1) M1E3(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 4
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E2(mask2); M3E2(mask3); M4E2(mask4)] ...
       [M2E3(mask2); M3E3(mask3); M4E3(mask4)] [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E1(mask1) M1E2(mask1) M1E3(mask1) M1E4(mask1)];
XGTtr = [M2GT(mask2); M3GT(mask3); M4GT(mask4)];








%%
%Train M1 M3 M4 Test M2
%Number of Energy 1
Xtr = [M1E1(mask1); M3E1(mask3); M4E1(mask4)];
Xte = M2E1(mask2);
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M1E2(mask1); M3E2(mask3); M4E2(mask4)];
Xte = M2E2(mask2);
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M1E3(mask1); M3E3(mask3); M4E3(mask4)];
Xte = M2E3(mask2);
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M1E4(mask1); M3E4(mask3); M4E4(mask4)];
Xte = M2E4(mask2);
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E2(mask1); M3E2(mask3); M4E2(mask4)]];
Xte = [M2E1(mask2) M2E2(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E3(mask1); M3E3(mask3); M4E3(mask4)]];
Xte = [M2E1(mask2) M2E3(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E1(mask2) M2E4(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E2(mask1); M3E2(mask3); M4E2(mask4)] [M1E3(mask1); M3E3(mask3); M4E3(mask4)]];
Xte = [M2E2(mask2) M2E3(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E2(mask1); M3E2(mask3); M4E2(mask4)] [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E2(mask2) M2E4(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E3(mask1); M3E3(mask3); M4E3(mask4)] [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E3(mask2) M2E4(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E2(mask1); M3E2(mask3); M4E2(mask4)]...
       [M1E3(mask1); M3E3(mask3); M4E3(mask4)]];
Xte = [M2E1(mask2) M2E2(mask2) M2E3(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E2(mask1); M3E2(mask3); M4E2(mask4)]...
       [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E1(mask2) M2E2(mask2) M2E3(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M1E2(mask1); M3E2(mask3); M4E2(mask4)] [M1E3(mask1); M3E3(mask3); M4E3(mask4)]...
       [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E2(mask2) M2E3(mask2) M2E4(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];
%%
%Number of Energy 4
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E2(mask1); M3E2(mask3); M4E2(mask4)]...
       [M1E3(mask1); M3E3(mask3); M4E3(mask4)] [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E1(mask2) M2E2(mask2) M2E3(mask2) M2E4(mask2)];
XGTtr = [M1GT(mask1); M3GT(mask3); M4GT(mask4)];




%%
%Train M1 M2 M4 Test M3
%Number of Energy 1
Xtr = [M1E1(mask1); M2E1(mask2); M4E1(mask4)];
Xte = M3E1(mask3);
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M1E2(mask1); M2E2(mask2); M4E2(mask4)];
Xte = M3E2(mask3);
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M1E3(mask1); M2E3(mask2); M4E3(mask4)];
Xte = M3E3(mask3);
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 1
Xtr = [M1E4(mask1); M2E4(mask2); M4E4(mask4)];
Xte = M3E4(mask3);
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E2(mask1); M2E2(mask2); M4E2(mask4)]];
Xte = [M3E1(mask3) M3E2(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E3(mask1); M2E3(mask2); M4E3(mask4)]];
Xte = [M3E1(mask3) M3E3(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E1(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E2(mask1); M2E2(mask2); M4E2(mask4)] [M1E3(mask1); M2E3(mask2); M4E3(mask4)]];
Xte = [M3E2(mask3) M3E3(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E2(mask1); M2E2(mask2); M4E2(mask4)] [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E2(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 2
Xtr = [[M1E3(mask1); M2E3(mask2); M4E3(mask4)] [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E3(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E2(mask1); M2E2(mask2); M4E2(mask4)]...
       [M1E3(mask1); M2E3(mask2); M4E3(mask4)]];
Xte = [M3E1(mask3) M3E2(mask3) M3E3(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E2(mask1); M2E2(mask2); M4E2(mask4)]...
       [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E1(mask3) M3E2(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 3
Xtr = [[M1E2(mask1); M2E2(mask2); M4E2(mask4)] [M1E3(mask1); M2E3(mask2); M4E3(mask4)]...
       [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E2(mask3) M3E3(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];
%%
%Number of Energy 4
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E2(mask1); M2E2(mask2); M4E2(mask4)] ...
       [M1E3(mask1); M2E3(mask2); M4E3(mask4)] [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E1(mask3) M3E2(mask3) M3E3(mask3) M3E4(mask3)];
XGTtr = [M1GT(mask1); M2GT(mask2); M4GT(mask4)];






%%
%Train M1 M2 M3 Test M4
%Number of Energy 1
Xtr = [M1E1(mask1); M2E1(mask2); M3E1(mask3)];
Xte = M4E1(mask4);
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 1
Xtr = [M1E2(mask1); M2E2(mask2); M3E2(mask3)];
Xte = M4E2(mask4);
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 1
Xtr = [M1E3(mask1); M2E3(mask2); M3E3(mask3)];
Xte = M4E3(mask4);
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 1
Xtr = [M1E4(mask1); M2E4(mask2); M3E4(mask3)];
Xte = M4E4(mask4);
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E2(mask1); M2E2(mask2); M3E2(mask3)]];
Xte = [M4E1(mask4) M4E2(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E3(mask1); M2E3(mask2); M3E3(mask3)]];
Xte = [M4E1(mask4) M4E3(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 2
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E1(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 2
Xtr = [[M1E2(mask1); M2E2(mask2); M3E2(mask3)] [M1E3(mask1); M2E3(mask2); M3E3(mask3)]];
Xte = [M4E2(mask4) M4E3(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 2
Xtr = [[M1E2(mask1); M2E2(mask2); M3E2(mask3)] [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E2(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 2
Xtr = [[M1E3(mask1); M2E3(mask2); M3E3(mask3)] [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E3(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 3
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E2(mask1); M2E2(mask2); M3E2(mask3)]...
       [M1E3(mask1); M2E3(mask2); M3E3(mask3)]];
Xte = [M4E1(mask4) M4E2(mask4) M4E3(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 3
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E2(mask1); M2E2(mask2); M3E2(mask3)]...
       [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E1(mask4) M4E2(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 3
Xtr = [[M1E2(mask1); M2E2(mask2); M3E2(mask3)] [M1E3(mask1); M2E3(mask2); M3E3(mask3)]...
       [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E2(mask4) M4E3(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];
%%
%Number of Energy 4
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E2(mask1); M2E2(mask2); M3E2(mask3)] ...
       [M1E3(mask1); M2E3(mask2); M3E3(mask3)] [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E1(mask4) M4E2(mask4) M4E3(mask4) M4E4(mask4)];
XGTtr = [M1GT(mask1); M2GT(mask2); M3GT(mask3)];





%%
%Number of Energy 1
for k = 1:K1
    tmp1 = Xtr(:,1);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.Sigma(:,:,k) = cov((tmp1(XGTtr == k)));
end
%%
%Number of Energy 2
for k = 1:K1
    tmp1 = Xtr(:,1); tmp2 = Xtr(:,2);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k)]));
end
%%
%Number of Energy 3
for k = 1:K1
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
%%
%Number of Energy 4
for k = 1:K1
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); tmp4 = Xtr(:,4);
    S.mu(k,1) = mean(tmp1(XGTtr == k));
    S.mu(k,2) = mean(tmp2(XGTtr == k));
    S.mu(k,3) = mean(tmp3(XGTtr == k));
    S.mu(k,4) = mean(tmp4(XGTtr == k));
    S.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k),tmp4(XGTtr == k)]));
end
%%
%atlas mouse1
atlas  = atlasfunc2(sig1,sig2,K,siz2,mask1,M2GT,M3GT);
%%
%atlas mouse2
atlas  = atlasfunc2(sig1,sig2,K,siz2,mask2,M1GT,M3GT);

%%
%atlas mouse3
atlas  = atlasfunc2(sig1,sig2,K,siz2,mask3,M1GT,M2GT);
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
I(mask3) = L;

%%
%Atlas_guided EM Mouse1
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,mask1,siz2);
JI= CalcuJI(Imap,M1GT,K-1);
disp("AtlasGuided EM_MAP result"); disp(JI);

%%
%Atlas_guided EM Mouse2
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,mask2,siz2);
JI= CalcuJI(Imap,M2GT,K-1);
disp("AtlasGuided EM_MAP result"); disp(JI);
%%
%Atlas_guided EM Mouse3
[Imap,L,PP,GMMMu,GMMSigma,GMMpro] = AtlasGuidedEM_kubo(Xte,atlas,S,K,mask3,siz2);
JI= CalcuJI(Imap,M3GT,K-1);
disp("AtlasGuided EM_MAP result"); disp(JI);