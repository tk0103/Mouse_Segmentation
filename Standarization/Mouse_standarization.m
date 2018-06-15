landmouse1 = [266,295,670; 307,281,670; 219,232,555; 296,209,555; 221,225,259; 245,211,256;];
landmouse2 = [298,260,644; 266,286,643; 349,304,526; 294,347,528; 327,299,236; 313,321,240;];
landmouse3 = [252,282,695; 295,272,696; 232,211,579; 300,191,580; 221,232,239; 246,220,236;];
%%
threold = 10;
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

IwM1E1 = apply_transformation_fast_3d( M1E1, invaffine1, siz );
IwM1E2 = apply_transformation_fast_3d( M1E2, invaffine1, siz );
IwM1E3 = apply_transformation_fast_3d( M1E3, invaffine1, siz );
IwM1E4 = apply_transformation_fast_3d( M1E4, invaffine1, siz );

IwM2E1 = apply_transformation_fast_3d( M2E1, invaffine2, siz );
IwM2E2 = apply_transformation_fast_3d( M2E2, invaffine2, siz );
IwM2E3 = apply_transformation_fast_3d( M2E3, invaffine2, siz );
IwM2E4 = apply_transformation_fast_3d( M2E4, invaffine2, siz );

IwM3E1 = apply_transformation_fast_3d( M3E1, invaffine3, siz );
IwM3E2 = apply_transformation_fast_3d( M3E2, invaffine3, siz );
IwM3E3 = apply_transformation_fast_3d( M3E3, invaffine3, siz );
IwM3E4 = apply_transformation_fast_3d( M3E4, invaffine3, siz );


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
pmask1 = mask1 - boneM1; pmask1 = logical(pmask1);
pmask2 = mask2 - boneM2; pmask2 = logical(pmask2);
pmask3 = mask3 - boneM3; pmask3 = logical(pmask3);

%%
%mask processing
pIwM1E1 = zeros(siz); pIwM1E2 = zeros(siz);  pIwM1E3 = zeros(siz);  pIwM1E4 = zeros(siz); pIwM1GT = zeros(siz);
pIwM1E1(pmask1) = IwM1E1(pmask1); pIwM1E2(pmask1) = IwM1E2(pmask1); pIwM1E3(pmask1) = IwM1E3(pmask1); pIwM1E4(pmask1) = IwM1E4(pmask1); pIwM1GT(pmask1) = IwM1GT(pmask1);

pIwM2E1 = zeros(siz); pIwM2E2 = zeros(siz);  pIwM2E3 = zeros(siz);  pIwM2E4 = zeros(siz); pIwM2GT = zeros(siz);
pIwM2E1(pmask2) = IwM2E1(pmask2); pIwM2E2(pmask2) = IwM2E2(pmask2); pIwM2E3(pmask2) = IwM2E3(pmask2); pIwM2E4(pmask2) = IwM2E4(pmask2); pIwM2GT(pmask2) = IwM2GT(pmask2);

pIwM3E1 = zeros(siz); pIwM3E2 = zeros(siz);  pIwM3E3 = zeros(siz);  pIwM3E4 = zeros(siz); pIwM3GT = zeros(siz);
pIwM3E1(pmask3) = IwM3E1(pmask3); pIwM3E2(pmask3) = IwM3E2(pmask3); pIwM3E3(pmask3) = IwM3E3(pmask3); pIwM3E4(pmask3) = IwM3E4(pmask3); pIwM3GT(pmask3) = IwM3GT(pmask3);

%%
imagesc(boneM1(:,:,325)');
axis tight equal
colormap gray
%%
%Save
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
