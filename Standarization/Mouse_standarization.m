%%
%Segmentation
M1E1 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\3DVolume\Mouse1\M1E1.raw','single');
M1E2 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\3DVolume\Mouse1\M1E2.raw','single');
M1E3 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\3DVolume\Mouse1\M1E3.raw','single');
M1E4 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\3DVolume\Mouse1\M1E4.raw','single');
M1GT = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\3DVolume\Mouse1\labelM1.raw','*uint8');
siz = [544,544,860];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz); 
M1E4 = reshape(M1E4,siz); M1GT = reshape(M1GT,siz);

%Input Mouse2
M2E1 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse2\M2E1.raw','single');
M2E2 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse2\M2E2.raw','single');
M2E3 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse2\M2E3.raw','single');
M2E4 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse2\M2E4.raw','single');
M2GT = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse2\labelM2.raw','*uint8');
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); 
M2E4 = reshape(M2E4,siz); M2GT = reshape(M2GT,siz);

%Input Mouse3
M3E1 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse3\M3E1.raw','single');
M3E2 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse3\M3E2.raw','single');
M3E3 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse3\M3E3.raw','single');
M3E4 = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse3\M3E4.raw','single');
M3GT = load_raw('\\tera\user\boku\NZ_data\Mouse\3DVolume\Mouse3\labelM3.raw','*uint8');
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); 
M3E4 = reshape(M3E4,siz); M3GT = reshape(M3GT,siz);

%%
%Input new data
M2E1 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\Wide_spectrum\Mouse2\E1-118.raw','*single');
M2E2 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\Wide_spectrum\Mouse2\E2-118.raw','*single');
%M2E3 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\Wide_spectrum\Mouse2\E3-118.raw','*single');
%M2E4 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\Wide_spectrum\Mouse2\E4-118.raw','*single');
siz = [672 672 792];
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz);
%M2E3 = reshape(M2E3,siz); 
%M2E4 = reshape(M2E4,siz);

M2E1 = M2E1(65:608,65:608,3:792);
M2E2 = M2E2(65:608,65:608,3:792);
M2E3 = M2E3(65:608,65:608,3:792);
M2E4 = M2E4(65:608,65:608,3:792);
%%
slice = 64;
subplot(1,2,1)
imagesc(M3E1(:,:,slice)');
axis tight equal 
colormap gray

subplot(1,2,2)
imagesc(M1E1(:,:,slice)');
axis tight equal 
colormap gray

%%
%registration
%sansho land3
%hudou land1 land2
%landmouse1 = [266,295,670; 307,281,670; 219,232,555; 296,209,555; 245,306,206; 337,266,206;];
%landmouse2 = [298,260,644; 266,286,643; 349,304,526; 294,347,528; 290,232,189; 234,309,207;];
%landmouse3 = [252,282,695; 295,272,696; 232,211,579; 300,191,580; 205,300,206; 301,253,191;];

%%
landmouse1 = [266,295,670; 307,281,670; 219,232,555; 296,209,555; 221,225,259; 245,211,256;];
landmouse2 = [298,260,644; 266,286,643; 349,304,526; 294,347,528; 327,299,236; 313,321,240;];
landmouse3 = [252,282,695; 295,272,696; 232,211,579; 300,191,580; 221,232,239; 246,220,236;];
%%
threold = 1000;
[~,Z1,transform1] = procrustes(landmouse3,landmouse1,'Reflection',false);
[~,Z2,transform2] = procrustes(landmouse3,landmouse2,'Reflection',false);
meanZ = (Z1 + Z2 +landmouse3 ) /3;

for round = 2:10
    disp(round);
    [~,Z1,transform1] = procrustes(meanZ,landmouse1,'Reflection',false);
    [~,Z2,transform2] = procrustes(meanZ,landmouse2,'Reflection',false);
    [~,Z3,transform3] = procrustes(meanZ,landmouse3,'Reflection',false);

    tmp = ((transform1.b + transform2.b + transform3.b)/3);
    transform1.b = transform1.b / tmp;
    transform2.b = transform2.b / tmp;
    transform3.b = transform3.b / tmp;
    
    Z1 = (transform1.b .* landmouse1 *transform1.T) + transform1.c;
    Z2 = (transform2.b .* landmouse2 *transform2.T) + transform2.c;
    Z3 = (transform3.b .* landmouse3 *transform3.T) + transform3.c;
    
    newmeanZ = (Z1+Z2+Z3)/3;
    thre = mean(sqrt(sum((meanZ - newmeanZ).^2,2)));
    
    if (abs(threold - thre)/threold) < 1e-3
        break;
    end

    threold = thre;
    meanZ = newmeanZ;
end
%%
affine1 = zeros(4,4);
affine1(1:3,1:3) = transform1.b .*transform1.T;
affine1(4,1:3) = transform1.c(1,:);
affine1(4,4) = 1;
affine1 = affine3d(affine1);
invaffine1 = invert(affine1);

affine2 = zeros(4,4);
affine2(1:3,1:3) = transform2.b .*transform2.T;
affine2(4,1:3) = transform2.c(1,:);
affine2(4,4) = 1;
affine2 = affine3d(affine2);
invaffine2 = invert(affine2);

affine3 = zeros(4,4);
affine3(1:3,1:3) = transform3.b .*transform3.T;
affine3(4,1:3) = transform3.c(1,:);
affine3(4,4) = 1;
affine3 = affine3d(affine3);
invaffine3 = invert(affine3);

%IwM1E1 = apply_transformation_fast_3d( M1E1, invaffine1, siz );
%IwM1E2 = apply_transformation_fast_3d( M1E2, invaffine1, siz );
%IwM1E3 = apply_transformation_fast_3d( M1E3, invaffine1, siz );
%IwM1E4 = apply_transformation_fast_3d( M1E4, invaffine1, siz );

IwM2E1 = apply_transformation_fast_3d( M2E1, invaffine2, siz );
IwM2E2 = apply_transformation_fast_3d( M2E2, invaffine2, siz );
%IwM2E3 = apply_transformation_fast_3d( M2E3, invaffine2, siz );
%IwM2E4 = apply_transformation_fast_3d( M2E4, invaffine2, siz );

%IwM3E1 = apply_transformation_fast_3d( M3E1, invaffine3, siz );
%IwM3E2 = apply_transformation_fast_3d( M3E2, invaffine3, siz );
%IwM3E3 = apply_transformation_fast_3d( M3E3, invaffine3, siz );
%IwM3E4 = apply_transformation_fast_3d( M3E4, invaffine3, siz );

%%
slice = 275;
subplot(1,3,1)
imagesc(test(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.7])

subplot(1,3,2)
imagesc(IwM3E1(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.7])

subplot(1,3,3)
imagesc(IwM1E1(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.7])

%%
%GT_transform
IwM1GT = zeros(siz); IwM2GT = zeros(siz); IwM3GT = zeros(siz);

for k =1:2
    mask = M1GT==k;
    nmask = not(mask);
    label = zeros(siz); nlabel = zeros(siz); 
    label(mask) = 1; nlabel(nmask) = 1;
    labeldist1 = bwdist(label); labeldist2 = bwdist(nlabel);
    
    labeldist1(nmask) = labeldist1(nmask)-0.5;
    labeldist2(mask) = -(labeldist2(mask)-0.5);
    labeldist = labeldist1 + labeldist2;
    
    tmp = apply_transformation_fast_3d( labeldist, invaffine1, siz );
    tmp = tmp<=0;
    IwM1GT(tmp) = k;
end
IwM1GT(IwM1GT==0) = 3;

for k =1:2
    mask = M2GT==k;
    nmask = not(mask);
    label = zeros(siz); nlabel = zeros(siz); 
    label(mask) = 1; nlabel(nmask) = 1;
    labeldist1 = bwdist(label); labeldist2 = bwdist(nlabel);
    
    labeldist1(nmask) = labeldist1(nmask)-0.5;
    labeldist2(mask) = -(labeldist2(mask)-0.5);
    labeldist = labeldist1 + labeldist2;
    
    tmp = apply_transformation_fast_3d( labeldist, invaffine2, siz );
    tmp = tmp<=0;
    IwM2GT(tmp) = k;
end
IwM2GT(IwM2GT==0) = 3;

for k =1:2
    mask = M3GT==k;
    nmask = not(mask);
    label = zeros(siz); nlabel = zeros(siz); 
    label(mask) = 1; nlabel(nmask) = 1;
    labeldist1 = bwdist(label); labeldist2 = bwdist(nlabel);
    
    labeldist1(nmask) = labeldist1(nmask)-0.5;
    labeldist2(mask) = -(labeldist2(mask)-0.5);
    labeldist = labeldist1 + labeldist2;
    
    tmp = apply_transformation_fast_3d( labeldist, invaffine3, siz );
    tmp = tmp<=0;
    IwM3GT(tmp) = k;
end
IwM3GT(IwM3GT==0) = 3;

%%
slice = 125;
subplot(1,2,1)
imagesc(pIwM1E1(:,:,slice)');
axis tight equal off
colormap default
%caxis([0 3])
colormap gray

subplot(1,2,2)
imagesc(test(:,:,slice)');
axis tight equal off
%caxis([0 3])

slice = 248;
subplot(1,3,3)
imagesc(pIwM2E1(:,:,slice)');
axis tight equal off
%caxis([0 3])


%%
save_raw(pIwM2E1,'C:\Users\yourb\Desktop\new2\pIwM1E1.raw','*double');
save_raw(pIwM2E2,'C:\Users\yourb\Desktop\new2\pIwM1E2.raw','*double');
%save_raw(te3,'C:\Users\yourb\Desktop\new2\pIwM1E3.raw','*double');
%save_raw(te4,'C:\Users\yourb\Desktop\new2\pIwM1E4.raw','*double');

%%
pmask2 = load_raw('\\tera\user\boku(�R�X�p���P��)\NZ_data\Mouse\Standarized_mouse_wide\pmaskM2.raw','*uint8');
pmask2 = reshape(pmask2,[544 544 860]);
pmask2 = pmask2(:,:,1:790);
%%
te1 = zeros(siz);
te1(pmask2) = pIwM2E1(pmask2);

te2 = zeros(siz);
te2(pmask2) = pIwM2E2(pmask2);

te3 = zeros(siz);
te3(pmask2) = pIwM2E3(pmask2);

te4 = zeros(siz);
te4(pmask2) = pIwM2E4(pmask2);
%%
imagesc(te1(:,:,332)');
%%
%background
M1bg = false(siz);
M1bg(IwM1E1 < graythresh(IwM1E1)) = 1;

M2bg = false(siz);
M2bg(IwM2E1 < graythresh(IwM2E1)) = 1;

M3bg = false(siz);
M3bg(IwM3E1 < graythresh(IwM3E1)) = 1;

%labeling & dilation
tube1 = false(siz);
L = bwlabeln(M1bg);
tube1(L==1) = true;
SE = strel('disk',16);
mask1 = imdilate(tube1,SE);
mask1 = or(mask1,M1bg);
mask1 = not(mask1);

tube2 = false(siz);
L = bwlabeln(M2bg);
tube2(L==1) = true;
mask2 = imdilate(tube2,SE);
mask2 = or(mask2,M2bg);
mask2 = not(mask2);

tube3 = false(siz);
L = bwlabeln(M3bg);
tube3(L==1) = true;
mask3 = imdilate(tube3,SE);
mask3 = or(mask3,M3bg);
mask3 = not(mask3);

%mask processing
pIwM1E1 = zeros(siz); pIwM1E2 = zeros(siz);  pIwM1E3 = zeros(siz);  pIwM1E4 = zeros(siz); pIwM1GT = zeros(siz);
pIwM1E1(mask1) = IwM1E1(mask1); pIwM1E2(mask1) = IwM1E2(mask1); pIwM1E3(mask1) = IwM1E3(mask1); pIwM1E4(mask1) = IwM1E4(mask1); pIwM1GT(mask1) = IwM1GT(mask1);

pIwM2E1 = zeros(siz); pIwM2E2 = zeros(siz);  pIwM2E3 = zeros(siz);  pIwM2E4 = zeros(siz); pIwM2GT = zeros(siz);
pIwM2E1(mask2) = IwM2E1(mask2); pIwM2E2(mask2) = IwM2E2(mask2); pIwM2E3(mask2) = IwM2E3(mask2); pIwM2E4(mask2) = IwM2E4(mask2); pIwM2GT(mask2) = IwM2GT(mask2);

pIwM3E1 = zeros(siz); pIwM3E2 = zeros(siz);  pIwM3E3 = zeros(siz);  pIwM3E4 = zeros(siz); pIwM3GT = zeros(siz);
pIwM3E1(mask3) = IwM3E1(mask3); pIwM3E2(mask3) = IwM3E2(mask3); pIwM3E3(mask3) = IwM3E3(mask3); pIwM3E4(mask3) = IwM3E4(mask3); pIwM3GT(mask3) = IwM3GT(mask3);

%%
nhood = zeros(3,3,3);
nhood(2,2,2) = 1.0;
for a = -1:1
    for b = -1:1
        for c= -1:1
            if(1<=abs(a)+abs(b)+abs(c) && abs(a)+abs(b)+abs(c)<=2)
                nhood(a+2,b+2,c+2) = 1.0;
            end
        end
    end
end
SE = strel('arbitrary',nhood);
%%
boneM1 = true(siz); boneM1ori = true(siz); tmp = true(siz);
boneM1ori(pIwM1E1 < graythresh(pIwM1E1(mask1))) = 0;
boneM1(pIwM1E1 < graythresh(pIwM1E1(mask1))) = 0;
diffM1 = pIwM1E2 - pIwM1E1;
tmp(diffM1 < 0.01) = 0;
boneM1 = boneM1 - tmp;
negative = find(boneM1<0);
boneM1(negative) = zeros(size(negative));

maxval = 4;
for n = 1:maxval
    boneM1 = imdilate(boneM1,SE);
end
for n = 1:maxval
    boneM1 = imerode(boneM1,SE);
end

BWM1 = zeros(siz);
L1 = bwconncomp(boneM1);
numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels);
BWM1(L1.PixelIdxList{idx}) = 1;
boneM1 = and(boneM1ori,BWM1);
%%
boneM2 = true(siz); boneM2ori = true(siz); tmp = true(siz);
boneM2ori(pIwM2E1 < graythresh(pIwM2E1(mask2))) = 0;
boneM2(pIwM2E1 < graythresh(pIwM2E1(mask2))) = 0;
diffM2 = pIwM2E2 - pIwM2E1;
tmp(diffM2 < 0.01) = 0;
boneM2 = boneM2 - tmp;
negative = find(boneM2<0);
boneM2(negative) = zeros(size(negative));

maxval = 2;
for n = 1:maxval
    boneM2 = imdilate(boneM2,SE);
end
for n = 1:maxval
    boneM2 = imerode(boneM2,SE);
end

BWM2 = zeros(siz);
L2 = bwconncomp(boneM2);
numPixels = cellfun(@numel,L2.PixelIdxList);
[~,idx] = max(numPixels);
BWM2(L2.PixelIdxList{idx}) = 1;
boneM2 = and(boneM2ori,BWM2);

%%
boneM3 = true(siz); boneM3ori = true(siz); tmp = true(siz);
boneM3ori(pIwM3E1 < graythresh(pIwM3E1(mask3))) = 0;
boneM3(pIwM3E1 < graythresh(pIwM3E1(mask3))) = 0;
diffM3 = pIwM3E2 - pIwM3E1;
tmp(diffM3 < 0.01) = 0;
boneM3 = boneM3 - tmp;
negative = find(boneM3<0);
boneM3(negative) = zeros(size(negative));

maxval = 4;
for n = 1:maxval
    boneM3 = imdilate(boneM3,SE);
end
for n = 1:maxval
    boneM3 = imerode(boneM3,SE);
end

BWM3 = zeros(siz);
L3 = bwconncomp(boneM3);
numPixels = cellfun(@numel,L3.PixelIdxList);
[~,idx] = max(numPixels);
BWM3(L3.PixelIdxList{idx}) = 1;
boneM3 = and(boneM3ori,BWM3);
%%
imagesc(pIwM1E1(:,:,220));
%%
pmask1 = mask1 - boneM1; pmask1 = logical(pmask1);
pmask2 = mask2 - boneM2; pmask2 = logical(pmask2);
pmask3 = mask3 - boneM3; pmask3 = logical(pmask3);
%%
imagesc(pIwM3E1(:,:,460)');
axis tight equal off
colormap gray
%caxis([0 0.5])
%%
%edge = [-0.05 -0.05:0.0005:0.1 0.1];
tmp =pIwM3E1(mask3);
histogram(tmp(:));
%%
%mask processing
pIwM1E1 = zeros(siz); pIwM1E2 = zeros(siz);  pIwM1E3 = zeros(siz);  pIwM1E4 = zeros(siz); pIwM1GT = zeros(siz);
pIwM1E1(pmask1) = IwM1E1(pmask1); pIwM1E2(pmask1) = IwM1E2(pmask1); pIwM1E3(pmask1) = IwM1E3(pmask1); pIwM1E4(pmask1) = IwM1E4(pmask1); pIwM1GT(pmask1) = IwM1GT(pmask1);

pIwM2E1 = zeros(siz); pIwM2E2 = zeros(siz);  pIwM2E3 = zeros(siz);  pIwM2E4 = zeros(siz); pIwM2GT = zeros(siz);
pIwM2E1(pmask2) = IwM2E1(pmask2); pIwM2E2(pmask2) = IwM2E2(pmask2); pIwM2E3(pmask2) = IwM2E3(pmask2); pIwM2E4(pmask2) = IwM2E4(pmask2); pIwM2GT(pmask2) = IwM2GT(pmask2);

pIwM3E1 = zeros(siz); pIwM3E2 = zeros(siz);  pIwM3E3 = zeros(siz);  pIwM3E4 = zeros(siz); pIwM3GT = zeros(siz);
pIwM3E1(pmask3) = IwM3E1(pmask3); pIwM3E2(pmask3) = IwM3E2(pmask3); pIwM3E3(pmask3) = IwM3E3(pmask3); pIwM3E4(pmask3) = IwM3E4(pmask3); pIwM3GT(pmask3) = IwM3GT(pmask3);

%%
imagesc(pIwM3E1(:,:,325)');
axis tight equal
colormap gray
%%
save_raw(pIwM1E1,'C:\Users\yourb\Desktop\new2\pIwM1E1.raw','*double');
save_raw(pIwM2E1,'C:\Users\yourb\Desktop\new2\pIwM2E1.raw','*double');
save_raw(pIwM3E1,'C:\Users\yourb\Desktop\new2\pIwM3E1.raw','*double');

save_raw(pIwM1E2,'C:\Users\yourb\Desktop\new2\pIwM1E2.raw','*double');
save_raw(pIwM2E2,'C:\Users\yourb\Desktop\new2\pIwM2E2.raw','*double');
save_raw(pIwM3E2,'C:\Users\yourb\Desktop\new2\pIwM3E2.raw','*double');

save_raw(pIwM1E3,'C:\Users\yourb\Desktop\new2\pIwM1E3.raw','*double');
save_raw(pIwM2E3,'C:\Users\yourb\Desktop\new2\pIwM2E3.raw','*double');
save_raw(pIwM3E3,'C:\Users\yourb\Desktop\new2\pIwM3E3.raw','*double');

save_raw(pIwM1E4,'C:\Users\yourb\Desktop\new2\pIwM1E4.raw','*double');
save_raw(pIwM2E4,'C:\Users\yourb\Desktop\new2\pIwM2E4.raw','*double');
save_raw(pIwM3E4,'C:\Users\yourb\Desktop\new2\pIwM3E4.raw','*double');

save_raw(pIwM1GT,'C:\Users\yourb\Desktop\new2\pIwM1GT.raw','*uint8');
save_raw(pIwM2GT,'C:\Users\yourb\Desktop\new2\pIwM2GT.raw','*uint8');
save_raw(pIwM3GT,'C:\Users\yourb\Desktop\new2\pIwM3GT.raw','*uint8');
%%
save_raw(pmask1,'C:\Users\yourb\Desktop\new2\pmaskM1.raw','*uint8');
save_raw(pmask2,'C:\Users\yourb\Desktop\new2\pmaskM2.raw','*uint8');
save_raw(pmask3,'C:\Users\yourb\Desktop\new2\pmaskM3.raw','*uint8');
%%
save_raw(boneM1,'C:\Users\yourb\Desktop\new2\boneM1.raw','*uint8');
save_raw(boneM2,'C:\Users\yourb\Desktop\new2\boneM2.raw','*uint8');
save_raw(boneM3,'C:\Users\yourb\Desktop\new2\boneM3.raw','*uint8');
%%
M3bgE1 = true(siz);
M3bgE1(pIwM3E1 < graythresh(pIwM3E1(mask3))) = 0;

M3bgE2 = true(siz);
M3bgE2(pIwM3E2 < graythresh(pIwM3E2(mask3))) = 0;

M3bgE3 = true(siz);
M3bgE3(pIwM3E3 < graythresh(pIwM3E3(mask3))) = 0;

M3bgE4 = true(siz);
M3bgE4(pIwM3E4 < graythresh(pIwM3E4(mask3))) = 0;
%%
aa = zeros(siz);
aa(M3bgbgE1) = IwM1E1(M3bgbgE1);
%%
imagesc(aa(:,:,260)');
colormap gray
%%
edge = [0.0 0.0:0.0001:0.6 0.6];
%histogram(pIwM3E1(M3bgE1),edge);
%histogram(pIwM3E2(M3bgE2),edge);
histogram(pM3E3(mask3),edge);
%histogram(pIwM3E4(M3bgE4),edge);
%%
figure;
histogram(pIwM3E1(M3bgE1),edge);
%%
M3bgbgE1 = true(siz);
M3bgbgE1(pIwM3E1 < 0.35) = 0;

%%
test = zeros(siz);
for c = 1:siz(3)
    for b = 1:siz(1)
        for a= 1:siz(2)
            if (M2GT(a,b,c)==100)
                test(a,b,c) = 1;
            end
            
            if (M2GT(a,b,c)==200)
                test(a,b,c) = 2;
            end
        end
    end
end
save_raw(test,'C:\Users\yourb\Desktop\new2\M2GT.raw','*uint8');

%%
%Segmentation
M1E1 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM1E1.raw','*double');
M1E2 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM1E2.raw','*double');
M1E3 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM1E3.raw','*double');
M1E4 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM1E4.raw','*double');
M1GT = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM1GT.raw','*uint8');
mask1 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pmaskM1.raw','*uint8');
siz = [544,544,860];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz); 
M1E4 = reshape(M1E4,siz); M1GT = reshape(M1GT,siz); mask1 = reshape(mask1,siz);
mask1 = logical(mask1);

%Input Mouse2
M2E1 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM2E1.raw','*double');
M2E2 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM2E2.raw','*double');
M2E3 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM2E3.raw','*double');
M2E4 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM2E4.raw','*double');
M2GT = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM2GT.raw','*uint8');
mask2 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pmaskM2.raw','*uint8');
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); 
M2E4 = reshape(M2E4,siz); M2GT = reshape(M2GT,siz); mask2 = reshape(mask2,siz);
mask2 = logical(mask2);

%Input Mouse3
M3E1 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM3E1.raw','*double');
M3E2 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM3E2.raw','*double');
M3E3 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM3E3.raw','*double');
M3E4 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM3E4.raw','*double');
M3GT = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pIwM3GT.raw','*uint8');
mask3 = load_raw('\\tera\user\boku\NZ_data\Mouse\Standarized_mouse\with_bone\pmaskM3.raw','*uint8');
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); 
M3E4 = reshape(M3E4,siz); M3GT = reshape(M3GT,siz); mask3 = reshape(mask3,siz);
mask3 = logical(mask3);
%%
imagesc(M2E1(:,:,210));
axis tight equal
%%
diff21M1 = M1E2 - M1E1;
diff21M2 = M2E2 - M2E1;
diff21M3 = M3E2 - M3E1;

SE = strel('sphere',3);
boneM1 = true(siz); tmp = true(siz); BWM1 = zeros(siz);
boneM1(M1E1 < 0.32) = 0; boneM1ori = boneM1;
tmp(diff21M1 < 0.015) = 0; boneM1(tmp) = 0;
boneM1 = imclose(boneM1,SE);
L1 = bwconncomp(boneM1);
numPixels = cellfun(@numel,L1.PixelIdxList);
[~,idx] = max(numPixels);
BWM1(L1.PixelIdxList{idx}) = 1;
boneM1 = and(boneM1ori,BWM1);
boneM1 = logical(boneM1);

%boneM2
SE = strel('sphere',3);
boneM2 = true(siz); tmp = true(siz); BWM2 = zeros(siz);
boneM2(M2E1 < 0.32) = 0; boneM2ori = boneM2;
tmp(diff21M2 < 0.015) = 0; boneM2(tmp) = 0;
boneM2 = imclose(boneM2,SE);
L2 = bwconncomp(boneM2);
numPixels = cellfun(@numel,L2.PixelIdxList);
[~,idx] = max(numPixels);
BWM2(L2.PixelIdxList{idx}) = 1;
boneM2 = and(boneM2ori,BWM2);
boneM2 = logical(boneM2);

%boneM3
SE = strel('sphere',5);
boneM3 = true(siz); tmp = true(siz); BWM3 = zeros(siz);
boneM3(M3E1 < 0.3) = 0; boneM3ori = boneM3;
tmp(diff21M3 < 0.015) = 0; boneM3(tmp) = 0;
boneM3 = imclose(boneM3,SE);
L3 = bwconncomp(boneM3);
numPixels = cellfun(@numel,L3.PixelIdxList);
[~,idx] = max(numPixels);
BWM3(L3.PixelIdxList{idx}) = 1;
boneM3 = and(boneM3ori,BWM3);
boneM3 = logical(boneM3);

boneM1 = not(boneM1);
boneM2 = not(boneM2);
boneM3 = not(boneM3);
%%
pM1E1 = zeros(siz); pM1E1(boneM1) = M1E1(boneM1); 
pM1E2 = zeros(siz); pM1E2(boneM1) = M1E2(boneM1); 
pM1E3 = zeros(siz); pM1E3(boneM1) = M1E3(boneM1); 
pM1E4 = zeros(siz); pM1E4(boneM1) = M1E4(boneM1); 
pM1GT = zeros(siz); pM1GT(boneM1) = M1GT(boneM1); 

pM2E1 = zeros(siz); pM2E1(boneM2) = M2E1(boneM2); 
pM2E2 = zeros(siz); pM2E2(boneM2) = M2E2(boneM2); 
pM2E3 = zeros(siz); pM2E3(boneM2) = M2E3(boneM2); 
pM2E4 = zeros(siz); pM2E4(boneM2) = M2E4(boneM2); 
pM2GT = zeros(siz); pM2GT(boneM2) = M2GT(boneM2); 

pM3E1 = zeros(siz); pM3E1(boneM3) = M3E1(boneM3); 
pM3E2 = zeros(siz); pM3E2(boneM3) = M3E2(boneM3); 
pM3E3 = zeros(siz); pM3E3(boneM3) = M3E3(boneM3); 
pM3E4 = zeros(siz); pM3E4(boneM3) = M3E4(boneM3); 
pM3GT = zeros(siz); pM3GT(boneM3) = M3GT(boneM3); 
%%
save_raw(pM1E1,'C:\Users\yourb\Desktop\new2\pIwM1E1.raw','*double');
save_raw(pM2E1,'C:\Users\yourb\Desktop\new2\pIwM2E1.raw','*double');
save_raw(pM3E1,'C:\Users\yourb\Desktop\new2\pIwM3E1.raw','*double');

save_raw(pM1E2,'C:\Users\yourb\Desktop\new2\pIwM1E2.raw','*double');
save_raw(pM2E2,'C:\Users\yourb\Desktop\new2\pIwM2E2.raw','*double');
save_raw(pM3E2,'C:\Users\yourb\Desktop\new2\pIwM3E2.raw','*double');

save_raw(pM1E3,'C:\Users\yourb\Desktop\new2\pIwM1E3.raw','*double');
save_raw(pM2E3,'C:\Users\yourb\Desktop\new2\pIwM2E3.raw','*double');
save_raw(pM3E3,'C:\Users\yourb\Desktop\new2\pIwM3E3.raw','*double');

save_raw(pM1E4,'C:\Users\yourb\Desktop\new2\pIwM1E4.raw','*double');
save_raw(pM2E4,'C:\Users\yourb\Desktop\new2\pIwM2E4.raw','*double');
save_raw(pM3E4,'C:\Users\yourb\Desktop\new2\pIwM3E4.raw','*double');

save_raw(pM1GT,'C:\Users\yourb\Desktop\new2\pIwM1GT.raw','*uint8');
save_raw(pM2GT,'C:\Users\yourb\Desktop\new2\pIwM2GT.raw','*uint8');
save_raw(pM3GT,'C:\Users\yourb\Desktop\new2\pIwM3GT.raw','*uint8');

save_raw(test1,'C:\Users\yourb\Desktop\new2\pmaskM1.raw','*uint8');
save_raw(test2,'C:\Users\yourb\Desktop\new2\pmaskM2.raw','*uint8');
save_raw(test3,'C:\Users\yourb\Desktop\new2\pmaskM3.raw','*uint8');

%%
boneM1 = not(boneM1);
boneM2 = not(boneM2);
boneM3 = not(boneM3);
%%
SE = strel('sphere',3);
testbone = imerode(boneM3,SE);
%%
imagesc(boneM3(:,:,295));
axis tight equal