%%
%train_mouse1 mouse2 test_mouse3
Xtr = [[pM1E2(pmask1); pM2E2(pmask2)] [pM1E3(pmask1); pM2E3(pmask2)] [pM1E4(pmask1); pM2E4(pmask2)] ];
Xte = [pM3E2(pmask3) pM3E3(pmask3) pM3E4(pmask3)];
XGTtr = [pM1GT(pmask1); pM2GT(pmask2)];
XGTte = pM3GT(pmask3);
%%
%initial_value
K=4;
sig1 = 5; %bladder
sig2 = 3; %kidneys

for k = 1:K
    tmp1 = Xtr(:,1);   tmp2 = Xtr(:,2);  tmp3 = Xtr(:,3); 
    SS.mu(k,1) = mean(tmp1(XGTtr == k));
    SS.mu(k,2) = mean(tmp2(XGTtr == k));
    SS.mu(k,3) = mean(tmp3(XGTtr == k));
    SS.Sigma(:,:,k) = cov(([tmp1(XGTtr == k),tmp2(XGTtr == k),tmp3(XGTtr == k)]));
end
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
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,lilelihood] = AtlasGuidedEM_kubo(Xte,atlas,SS,K,pmask3,siz2,30);
JI= CalcuJI(Imap,pM3GT,K-1);
disp("EM_MAP result")
disp(JI);
%%
temp = Imap == 1; radi = power(bwarea(temp(:))*6/pi,1/3);
D = bwdist(not(temp),'euclidean'); 
[~,I] =  max(D(:)); [x,y,z] = ind2sub(siz2,I); coor(1,:) = [x y z radi]; 

temp = Imap == 2; radi = power(bwarea(temp(:))*6/pi,1/3);
D = bwdist(not(temp),'euclidean'); 
[~,I] =  max(D(:)); [x,y,z] = ind2sub(siz2,I); coor(2,:) = [x y z radi]; 

temp = Imap == 3; radi = power(bwarea(temp(:))*6/pi,1/3);
D = bwdist(not(temp),'euclidean'); 
[~,I] =  max(D(:)); [x,y,z] = ind2sub(siz2,I); coor(3,:) = [x y z radi]; 
clearvars x y z radi I D temp
%%
newmaskM3 = zeros(siz2);
new