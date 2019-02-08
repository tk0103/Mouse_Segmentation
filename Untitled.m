In = pIwM2E2; th = 0.32;
siz = size(In);
tmp1 = pIwM1GT; tmp1(tmp1 == 4) = 0; tmp1 = logical(tmp1);
tmp2 = pIwM3GT; tmp2(tmp2 == 4) = 0; tmp2 = logical(tmp2);
tmp3 = pIwM4GT; tmp3(tmp3 == 4) = 0; tmp3 = logical(tmp3);
GTtemp = (tmp1 + tmp2 + tmp3)/3.0;
SE = strel('sphere',2); 
%%
boneth = true(siz);
boneth(In < th) = 0;
boneth  =imopen(boneth,SE);
tmp = zeros(siz); tmp(In == 0) = 1;
L = bwconncomp(tmp);
mask = zeros(siz);
[~,idx] = max(cellfun(@numel,L.PixelIdxList));
mask(L.PixelIdxList{idx}) = 1;
distIm = bwdist(logical(mask));
CC = bwconncomp(boneth,26);
L = labelmatrix(CC);
%%
dist = zeros(1,CC.NumObjects);
atval = zeros(1,CC.NumObjects);
stats = regionprops(CC,'Centroid','Area');
for n = 1:CC.NumObjects
    tmp = zeros(siz);
    tmp(L == n)  = 1;
    tmp = logical(tmp);
    dist(n) = mean(distIm(tmp));
    atval(n) = mean(GTtemp(tmp));
end
%%
out = ismember(L, find( atval<=0.34  &...
    dist<=50 & [stats.Area]>= median([stats.Area])));
out = imdilate(out,SE);

%%
a = [stats.Area];
scatter3(a,dist,atval,'MarkerEdgeColor','none','MarkerFaceColor','r');
xlabel('体積（ボクセル数）')
ylabel('体表からの距離')
zlabel('確率値');
%%
map2 = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.4
    0.8, 0.7, 0.4
    0.6 , 0.5, 0.5
    0.9, 0.9, 0];
imagesc(L(:,:,70)');
caxis([0 6])
colormap(map2)
axis tight equal off
%caxis([0 0.7])
%%
imagesc(tmp1(:,:,280)');
axis tight equal
%caxis([0 0.7])
%%
save_raw(out,'C:\Users\yourb\Desktop\new3\boneM4.raw','*uint8');

%%
boneM1 = load_raw('C:\Users\yourb\Desktop\new3\boneM1.raw','*uint8');
boneM2 = load_raw('C:\Users\yourb\Desktop\new3\boneM2.raw','*uint8');
boneM3 = load_raw('C:\Users\yourb\Desktop\new3\boneM3.raw','*uint8');
boneM4 = load_raw('C:\Users\yourb\Desktop\new3\boneM4.raw','*uint8');
boneM1 = reshape(boneM1,siz);
boneM2 = reshape(boneM2,siz);
boneM3 = reshape(boneM3,siz);
boneM4 = reshape(boneM4,siz);

%%
imagesc(boneM3(:,:,200)');
axis tight equal off
%colormap(gray)
%caxis([0 0.7])
%%
imagesc(M1E1(:,:,70)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
pmask1 = mask1 - logical(boneM1);
pmask2 = mask2 - logical(boneM2);
pmask3 = mask3 - logical(boneM3);
pmask4 = mask4 - logical(boneM4);

pmask1 = logical(pmask1);
pmask2 = logical(pmask2);
pmask3 = logical(pmask3);
pmask4 = logical(pmask4);
%%
pM1E1 = zeros(siz);
pM1E1(pmask1) = pIwM1E1(pmask1);
pM1E2 = zeros(siz);
pM1E2(pmask1) = pIwM1E2(pmask1);
pM1E3 = zeros(siz);
pM1E3(pmask1) = pIwM1E3(pmask1);
pM1E4 = zeros(siz);
pM1E4(pmask1) = pIwM1E4(pmask1);
% pwM1 = zeros(siz);
% pwM1(pmask1) = pIwM1(pmask1);

pM2E1 = zeros(siz);
pM2E1(pmask2) = pIwM2E1(pmask2);
pM2E2 = zeros(siz);
pM2E2(pmask2) = pIwM2E2(pmask2);
pM2E3 = zeros(siz);
pM2E3(pmask2) = pIwM2E3(pmask2);
pM2E4 = zeros(siz);
pM2E4(pmask2) = pIwM2E4(pmask2);
% pwM2 = zeros(siz);
% pwM2(pmask2) = pIwM2(pmask2);

pM3E1 = zeros(siz);
pM3E1(pmask3) = pIwM3E1(pmask3);
pM3E2 = zeros(siz);
pM3E2(pmask3) = pIwM3E2(pmask3);
pM3E3 = zeros(siz);
pM3E3(pmask3) = pIwM3E3(pmask3);
pM3E4 = zeros(siz);
pM3E4(pmask3) = pIwM3E4(pmask3);
% pwM3 = zeros(siz);
% pwM3(pmask3) = pIwM3(pmask3);

pM4E1 = zeros(siz);
pM4E1(pmask4) = pIwM4E1(pmask4);
pM4E2 = zeros(siz);
pM4E2(pmask4) = pIwM4E2(pmask4);
pM4E3 = zeros(siz);
pM4E3(pmask4) = pIwM4E3(pmask4);
pM4E4 = zeros(siz);
pM4E4(pmask4) = pIwM4E4(pmask4);
% pwM4 = zeros(siz);
% pwM4(pmask4) = pIwM4(pmask4);


pM1GT = zeros(siz);
pM1GT(pmask1) = pIwM1GT(pmask1);
pM2GT = zeros(siz);
pM2GT(pmask2) = pIwM2GT(pmask2);
pM3GT = zeros(siz);
pM3GT(pmask3) = pIwM3GT(pmask3);
pM4GT = zeros(siz);
pM4GT(pmask4) = pIwM4GT(pmask4);
%%
M1E1 = pM1E1;
M1E2 = pM1E2;
M1E3 = pM1E3;
M1E4 = pM1E4;

M2E1 = pM2E1;
M2E2 = pM2E2;
M2E3 = pM2E3;
M2E4 = pM2E4;

M3E1 = pM3E1;
M3E2 = pM3E2;
M3E3 = pM3E3;
M3E4 = pM3E4;

M4E1 = pM4E1;
M4E2 = pM4E2;
M4E3 = pM4E3;
M4E4 = pM4E4;

M1GT = pM1GT;
M2GT = pM2GT;
M3GT = pM3GT;
M4GT = pM4GT;

mask1 = pmask1;
mask2 = pmask2;
mask3 = pmask3;
mask4 = pmask4;


wM1E1 = pwM1;
wM2E1 = pwM2;
wM3E1 = pwM3;
wM4E1 = pwM4;
%%
temp = zeros(siz);
temp(pM1GT == 1) = 8;
SE = strel('sphere',1); 
%%
temp2 = imerode(temp,SE);
%%
pM1GT(pM1GT == 1) = 4;
pM1GT(temp2 == 8) = 1;
%%
temp2 = temp- temp2; 
%%
imagesc(pIwM3E2(:,:,70)');
axis tight equal off
colormap(gray)
caxis([0 0.7])
%%
imagesc(pmask3(:,:,70)');
axis tight equal off
%colormap(map)
%caxis([0 4])
%%
pM1GT = pM1GT - uint8(temp2);
%%
M1GT(M1GT == 5) = 1;
M1GT = uint8(M1GT);