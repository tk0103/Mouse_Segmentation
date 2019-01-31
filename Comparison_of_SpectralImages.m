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
Xtr = [[M2E1(mask2); M3E1(mask3); M4E1(mask4)] [M2E3(mask2); M3E3(mask3);M4E3(mask4)] ...
      [M2E4(mask2); M3E4(mask3); M4E4(mask4)]];
Xte = [M1E1(mask1) M1E3(mask1) M1E4(mask1)];
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
Xtr = [[M1E1(mask1); M3E1(mask3); M4E1(mask4)] [M1E3(mask1); M3E3(mask3); M4E3(mask4)]...
       [M1E4(mask1); M3E4(mask3); M4E4(mask4)]];
Xte = [M2E1(mask2) M2E3(mask2) M2E4(mask2)];
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
Xtr = [[M1E1(mask1); M2E1(mask2); M4E1(mask4)] [M1E3(mask1); M2E3(mask2); M4E3(mask4)]...
       [M1E4(mask1); M2E4(mask2); M4E4(mask4)]];
Xte = [M3E1(mask3) M3E3(mask3) M3E4(mask3)];
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
imagesc(M4E2(:,:,100)');
axis tight equal




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
Xtr = [[M1E1(mask1); M2E1(mask2); M3E1(mask3)] [M1E3(mask1); M2E3(mask2); M3E3(mask3)]...
       [M1E4(mask1); M2E4(mask2); M3E4(mask3)]];
Xte = [M4E1(mask4) M4E3(mask4) M4E4(mask4)];
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
%Mouse1
%Number of Energy 1
Xtebla  = M1E1(blamask);
XteLkid = M1E1(Lkidmask);
XteRkid = M1E1(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M1E2(blamask);
XteLkid = M1E2(Lkidmask);
XteRkid = M1E2(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M1E3(blamask);
XteLkid = M1E3(Lkidmask);
XteRkid = M1E3(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M1E4(blamask);
XteLkid = M1E4(Lkidmask);
XteRkid = M1E4(Rkidmask);
%%
%Number of Energy 2
Xtebla  = [M1E1(blamask)  M1E2(blamask)];
XteLkid = [M1E1(Lkidmask) M1E2(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E2(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M1E1(blamask)  M1E3(blamask)];
XteLkid = [M1E1(Lkidmask) M1E3(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M1E1(blamask)  M1E4(blamask)];
XteLkid = [M1E1(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M1E2(blamask)  M1E3(blamask)];
XteLkid = [M1E2(Lkidmask) M1E3(Lkidmask)];
XteRkid = [M1E2(Rkidmask) M1E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M1E2(blamask)  M1E4(blamask)];
XteLkid = [M1E2(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E2(Rkidmask) M1E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M1E3(blamask)  M1E4(blamask)];
XteLkid = [M1E3(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E3(Rkidmask) M1E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M1E1(blamask)  M1E2(blamask)  M1E3(blamask)];
XteLkid = [M1E1(Lkidmask) M1E2(Lkidmask) M1E3(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E2(Rkidmask) M1E3(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M1E1(blamask)  M1E2(blamask)  M1E4(blamask)];
XteLkid = [M1E1(Lkidmask) M1E2(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E2(Rkidmask) M1E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M1E2(blamask)  M1E3(blamask)  M1E4(blamask)];
XteLkid = [M1E2(Lkidmask) M1E3(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E2(Rkidmask) M1E3(Rkidmask) M1E4(Rkidmask)];
%%
%Number of Energy 4
Xtebla  = [M1E1(blamask)  M1E2(blamask)  M1E3(blamask)  M1E4(blamask)];
XteLkid = [M1E1(Lkidmask) M1E2(Lkidmask) M1E3(Lkidmask) M1E4(Lkidmask)];
XteRkid = [M1E1(Rkidmask) M1E2(Rkidmask) M1E3(Rkidmask) M1E4(Rkidmask)];





%%
%Mouse2
%Number of Energy 1
Xtebla  = M2E1(blamask);
XteLkid = M2E1(Lkidmask);
XteRkid = M2E1(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M2E2(blamask);
XteLkid = M2E2(Lkidmask);
XteRkid = M2E2(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M2E3(blamask);
XteLkid = M2E3(Lkidmask);
XteRkid = M2E3(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M2E4(blamask);
XteLkid = M2E4(Lkidmask);
XteRkid = M2E4(Rkidmask);
%%
%Number of Energy 2
Xtebla  = [M2E1(blamask)  M2E2(blamask)];
XteLkid = [M2E1(Lkidmask) M2E2(Lkidmask)];
XteRkid = [M2E1(Rkidmask) M2E2(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M2E1(blamask)  M2E3(blamask)];
XteLkid = [M2E1(Lkidmask) M2E3(Lkidmask)];
XteRkid = [M2E1(Rkidmask) M2E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M2E1(blamask)  M2E4(blamask)];
XteLkid = [M2E1(Lkidmask) M2E4(Lkidmask)];
XteRkid = [M2E1(Rkidmask) M2E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M2E2(blamask)  M2E3(blamask)];
XteLkid = [M2E2(Lkidmask) M2E3(Lkidmask)];
XteRkid = [M2E2(Rkidmask) M2E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M2E2(blamask)  M2E4(blamask)];
XteLkid = [M2E2(Lkidmask) M2E4(Lkidmask)];
XteRkid = [M2E2(Rkidmask) M2E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M2E3(blamask)  M2E4(blamask)];
XteLkid = [M2E3(Lkidmask) M2E4(Lkidmask)];
XteRkid = [M2E3(Rkidmask) M2E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M2E1(blamask)  M2E2(blamask)  M2E3(blamask)];
XteLkid = [M2E1(Lkidmask) M2E2(Lkidmask) M2E3(Lkidmask)];
XteRkid = [M2E1(Rkidmask) M2E2(Rkidmask) M2E3(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M2E1(blamask)  M2E2(blamask)  M2E4(blamask)];
XteLkid = [M2E1(Lkidmask) M2E2(Lkidmask) M2E4(Lkidmask)];
XteRkid = [M2E1(Rkidmask) M2E2(Rkidmask) M2E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M2E2(blamask)  M2E3(blamask)  M2E4(blamask)];
XteLkid = [M2E2(Lkidmask) M2E3(Lkidmask) M2E4(Lkidmask)];
XteRkid = [M2E2(Rkidmask) M2E3(Rkidmask) M2E4(Rkidmask)];
%%
%Number of Energy 4
Xtebla  = [M2E1(blamask)  ME2(blamask)   M2E3(blamask)  M2E4(blamask)];
XteLkid = [M2E1(Lkidmask) M2E2(Lkidmask) M2E3(Lkidmask) M2E4(Lkidmask)];
XteRkid = [M2E1(Rkidmask) M2E2(Rkidmask) M2E3(Rkidmask) M2E4(Rkidmask)];





%%
%Mouse3
%Number of Energy 1
Xtebla  = M3E1(blamask);
XteLkid = M3E1(Lkidmask);
XteRkid = M3E1(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M3E2(blamask);
XteLkid = M3E2(Lkidmask);
XteRkid = M3E2(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M3E3(blamask);
XteLkid = M3E3(Lkidmask);
XteRkid = M3E3(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M3E4(blamask);
XteLkid = M3E4(Lkidmask);
XteRkid = M3E4(Rkidmask);
%%
%Number of Energy 2
Xtebla  = [M3E1(blamask)  M3E2(blamask)];
XteLkid = [M3E1(Lkidmask) M3E2(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E2(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M3E1(blamask)  M3E3(blamask)];
XteLkid = [M3E1(Lkidmask) M3E3(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M3E1(blamask)  M3E4(blamask)];
XteLkid = [M3E1(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M3E2(blamask)  M3E3(blamask)];
XteLkid = [M3E2(Lkidmask) M3E3(Lkidmask)];
XteRkid = [M3E2(Rkidmask) M3E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M3E2(blamask)  M3E4(blamask)];
XteLkid = [M3E2(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E2(Rkidmask) M3E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M3E3(blamask)  M3E4(blamask)];
XteLkid = [M3E3(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E3(Rkidmask) M3E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M3E1(blamask)  M3E2(blamask)  M3E3(blamask)];
XteLkid = [M3E1(Lkidmask) M3E2(Lkidmask) M3E3(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E2(Rkidmask) M3E3(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M3E1(blamask)  M3E2(blamask)  M3E4(blamask)];
XteLkid = [M3E1(Lkidmask) M3E2(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E2(Rkidmask) M3E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M3E2(blamask)  M3E3(blamask)  M3E4(blamask)];
XteLkid = [M3E2(Lkidmask) M3E3(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E2(Rkidmask) M3E3(Rkidmask) M3E4(Rkidmask)];
%%
%Number of Energy 4
Xtebla  = [M3E1(blamask)  M3E2(blamask)  M3E3(blamask)  M3E4(blamask)];
XteLkid = [M3E1(Lkidmask) M3E2(Lkidmask) M3E3(Lkidmask) M3E4(Lkidmask)];
XteRkid = [M3E1(Rkidmask) M3E2(Rkidmask) M3E3(Rkidmask) M3E4(Rkidmask)];





%%
%Mouse4
%Number of Energy 1
Xtebla  = M4E1(blamask);
XteLkid = M4E1(Lkidmask);
XteRkid = M4E1(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M4E2(blamask);
XteLkid = M4E2(Lkidmask);
XteRkid = M4E2(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M4E3(blamask);
XteLkid = M4E3(Lkidmask);
XteRkid = M4E3(Rkidmask);
%%
%Number of Energy 1
Xtebla  = M4E4(blamask);
XteLkid = M4E4(Lkidmask);
XteRkid = M4E4(Rkidmask);
%%
%Number of Energy 2
Xtebla  = [M4E1(blamask)  M4E2(blamask)];
XteLkid = [M4E1(Lkidmask) M4E2(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E2(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M4E1(blamask)  M4E3(blamask)];
XteLkid = [M4E1(Lkidmask) M4E3(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M4E1(blamask)  M4E4(blamask)];
XteLkid = [M4E1(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M4E2(blamask)  M4E3(blamask)];
XteLkid = [M4E2(Lkidmask) M4E3(Lkidmask)];
XteRkid = [M4E2(Rkidmask) M4E3(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M4E2(blamask)  M4E4(blamask)];
XteLkid = [M4E2(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E2(Rkidmask) M4E4(Rkidmask)];
%%
%Number of Energy 2
Xtebla  = [M4E3(blamask)  M4E4(blamask)];
XteLkid = [M4E3(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E3(Rkidmask) M4E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M4E1(blamask)  M4E2(blamask)  M4E3(blamask)];
XteLkid = [M4E1(Lkidmask) M4E2(Lkidmask) M4E3(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E2(Rkidmask) M4E3(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M4E1(blamask)  M4E2(blamask)  M4E4(blamask)];
XteLkid = [M4E1(Lkidmask) M4E2(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E2(Rkidmask) M4E4(Rkidmask)];
%%
%Number of Energy 3
Xtebla  = [M4E2(blamask)  M4E3(blamask)  M4E4(blamask)];
XteLkid = [M4E2(Lkidmask) M4E3(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E2(Rkidmask) M4E3(Rkidmask) M4E4(Rkidmask)];
%%
%Number of Energy 4
Xtebla  = [M4E1(blamask)  M4E2(blamask)  M4E3(blamask)  M4E4(blamask)];
XteLkid = [M4E1(Lkidmask) M4E2(Lkidmask) M4E3(Lkidmask) M4E4(Lkidmask)];
XteRkid = [M4E1(Rkidmask) M4E2(Rkidmask) M4E3(Rkidmask) M4E4(Rkidmask)];