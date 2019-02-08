x = 0; y = 0;
landmouse1 = [257,257,780;  307,281,670; 266,295,670; 296,209,555; 219,232,555; 221,225,259; 245,211,256;];
landmouse2 = [257,257,753;  266,286,643; 298,260,644; 294,347,535; 349,304,531; 327,299,266; 313,321,240;];
landmouse3 = [257,257,808;  295,272,696; 252,282,695; 300,191,580; 232,211,579; 221,232,239; 246,220,236;];
landmouse4 = [257,257,827;  291,257,700; 257,284,699;  271+x,192+y,591; 210+x,215+y,584; 212,273,205; 301,181,205;];
%landmouse4 = [ 260,290,699;285,257,700; 200,261,584; 231,190,591;];
%%
x = 0; y = 0;
landmouse1 = [257,257,780;  307,281,670; 266,295,670; 296,209,555; 219,232,555; 221,225,259; 245,211,256; ];
landmouse2 = [257,257,733;  266,286,643; 298,260,644; 294,347,535; 349,304,531; 327,299,236; 313,321,210;];
landmouse3 = [257,257,808;  295,272,696; 252,282,695; 300,191,580; 232,211,579; 221,232,239; 246,220,236;];
landmouse4 = [257,257,827;  291,257,700; 257,284,699;  271+x,192+y,591; 210+x,215+y,584;212,273,205; 301,181,205;];
%landmouse4 = [ 260,290,699;285,257,700; 200,261,584; 231,190,591;];
%%
threold = 10;
[~,Z1,transform1] = procrustes(landmouse3,landmouse1,'Reflection',false);
[~,Z2,transform2] = procrustes(landmouse3,landmouse2,'Reflection',false);
[~,Z3,transform3] = procrustes(landmouse3,landmouse4,'Reflection',false);
meanZ = (Z1 + Z2 + Z3 + landmouse3 ) /4;

for round = 2:10
    disp(round);
    [~,Z1,transform1] = procrustes(meanZ,landmouse1,'Reflection',false);
    [~,Z2,transform2] = procrustes(meanZ,landmouse2,'Reflection',false);
    [~,Z3,transform3] = procrustes(meanZ,landmouse3,'Reflection',false);
    [~,Z4,transform4] = procrustes(meanZ,landmouse4,'Reflection',false);
    
    tmp = ((transform1.b + transform2.b + transform3.b + transform4.b)/4);
    transform1.b = transform1.b / tmp;
    transform2.b = transform2.b / tmp;
    transform3.b = transform3.b / tmp;
    transform4.b = transform4.b / tmp;
    
    Z1 = (transform1.b .* landmouse1 *transform1.T) + transform1.c;
    Z2 = (transform2.b .* landmouse2 *transform2.T) + transform2.c;
    Z3 = (transform3.b .* landmouse3 *transform3.T) + transform3.c;
    Z4 = (transform4.b .* landmouse4 *transform4.T) + transform4.c;
    
    newmeanZ = (Z1 + Z2 + Z3 + Z4) / 4;
    thre = mean(sqrt(sum((meanZ - newmeanZ).^2,2)));
    
    if (abs(threold - thre)/threold) < 1e-3
        break;
    end

    threold = thre;
    meanZ = newmeanZ;
end

affine1 = zeros(4,4);
affine1(1:3,1:3) = transform1.b .*transform1.T;
affine1(4,1:3) = transform1.c(1,:);
affine1(4,4) = 1;
invaffine1 = invert(affine3d(affine1));

affine2 = zeros(4,4);
affine2(1:3,1:3) = transform2.b .*transform2.T;
affine2(4,1:3) = transform2.c(1,:);
affine2(4,4) = 1;
invaffine2 = invert(affine3d(affine2));

affine3 = zeros(4,4);
affine3(1:3,1:3) = transform3.b .*transform3.T;
affine3(4,1:3) = transform3.c(1,:);
affine3(4,4) = 1;
invaffine3 = invert(affine3d(affine3));

affine4 = zeros(4,4);
affine4(1:3,1:3) = transform4.b .*transform4.T;
affine4(4,1:3) = transform4.c(1,:);
affine4(4,4) = 1;
invaffine4 = invert(affine3d(affine4));
%%
IwM1E1 = apply_transformation_fast_3d( M1E1, invaffine1, siz );
IwM1E2 = apply_transformation_fast_3d( M1E2, invaffine1, siz );
IwM1E3 = apply_transformation_fast_3d( M1E3, invaffine1, siz );
IwM1E4 = apply_transformation_fast_3d( M1E4, invaffine1, siz );
%IwM1 = apply_transformation_fast_3d( wM1, invaffine1, siz );

IwM2E1 = apply_transformation_fast_3d( M2E1, invaffine2, siz );
IwM2E2 = apply_transformation_fast_3d( M2E2, invaffine2, siz );
IwM2E3 = apply_transformation_fast_3d( M2E3, invaffine2, siz );
IwM2E4 = apply_transformation_fast_3d( M2E4, invaffine2, siz );
%IwM2 = apply_transformation_fast_3d( wM2, invaffine2, siz );

IwM3E1 = apply_transformation_fast_3d( M3E1, invaffine3, siz );
IwM3E2 = apply_transformation_fast_3d( M3E2, invaffine3, siz );
IwM3E3 = apply_transformation_fast_3d( M3E3, invaffine3, siz );
IwM3E4 = apply_transformation_fast_3d( M3E4, invaffine3, siz );
%IwM3 = apply_transformation_fast_3d( wM3, invaffine3, siz );

IwM4E1 = apply_transformation_fast_3d( M4E1, invaffine4, siz );
IwM4E2 = apply_transformation_fast_3d( M4E2, invaffine4, siz );
IwM4E3 = apply_transformation_fast_3d( M4E3, invaffine4, siz );
IwM4E4 = apply_transformation_fast_3d( M4E4, invaffine4, siz );
%IwM4 = apply_transformation_fast_3d( wM4, invaffine4, siz );
%%
slice = 419;
subplot(2,2,1)
imagesc(IwM1E2(:,:,slice)');
axis tight equal off
colormap(gray)
caxis([0 0.7])

subplot(2,2,2)
imagesc(IwM2E2(:,:,slice)'); 
axis tight equal off
colormap(gray)
caxis([0 0.7])

subplot(2,2,3)
imagesc(IwM3E2(:,:,slice)');
axis tight equal off
colormap(gray)
caxis([0 0.7])

subplot(2,2,4)
imagesc(IwM4E2(:,:,slice)');
axis tight equal off
colormap(gray)
caxis([0 0.7])
%%
save_raw(IwM4E2,'C:\Users\yourb\Desktop\pIwM4E2.raw','*single');
%%
tm =(IwM4E2 + IwM2E2 + IwM3E2 +IwM1E2)/4; 
%%
slice = 250;
imagesc(mask1(:,:,slice)');
axis tight equal off
colormap(gray)
caxis([0 0.7])
%%
%GT_transform
IwM1GT = zeros(siz);
for k =1:3
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
IwM1GT(IwM1GT==0) = 4;

IwM2GT = zeros(siz);
for k =1:3
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
IwM2GT(IwM2GT==0) = 4;

IwM3GT = zeros(siz);
for k =1:3
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
IwM3GT(IwM3GT==0) = 4;

IwM4GT = zeros(siz);
for k =1:3
    mask = M4GT==k;
    nmask = not(mask);
    label = zeros(siz); nlabel = zeros(siz); 
    label(mask) = 1; nlabel(nmask) = 1;
    labeldist1 = bwdist(label); labeldist2 = bwdist(nlabel);
    
    labeldist1(nmask) = labeldist1(nmask)-0.5;
    labeldist2(mask) = -(labeldist2(mask)-0.5);
    labeldist = labeldist1 + labeldist2;
    
    tmp = apply_transformation_fast_3d( labeldist, invaffine4, siz );
    tmp = tmp<=0;
    IwM4GT(tmp) = k;
end
IwM4GT(IwM4GT==0) = 4;

%%
ind = find(IwM1GT == 1); [~,~,z] = ind2sub(siz,ind);
minmunZ(1) = min(z);

ind = find(IwM2GT == 1); [~,~,z] = ind2sub(siz,ind);
minmunZ(2) = min(z);

ind = find(IwM3GT == 1); [~,~,z] = ind2sub(siz,ind);
minmunZ(3) = min(z);

ind = find(IwM4GT == 1); [~,~,z] = ind2sub(siz,ind);
minmunZ(4) = min(z);
st = min(minmunZ);
disp(min(minmunZ));

ind = find(IwM1GT == 3); [~,~,z] = ind2sub(siz,ind);
maximumZ(1) = max(z);

ind = find(IwM2GT == 3); [~,~,z] = ind2sub(siz,ind);
maximumZ(2) = max(z);

ind = find(IwM3GT == 3); [~,~,z] = ind2sub(siz,ind);
maximumZ(3) = max(z);

ind = find(IwM4GT == 3); [~,~,z] = ind2sub(siz,ind);
maximumZ(4) = max(z);
disp(max(maximumZ));
en  = max(maximumZ);
%%
%background
M1bg = false(siz);
M1bg(IwM1E1 < graythresh(IwM1E1)) = 1;

M2bg = false(siz);
M2bg(IwM2E1 < graythresh(IwM2E1)) = 1;

M3bg = false(siz);
M3bg(IwM3E1 < graythresh(IwM3E1)) = 1;

M4bg = false(siz);
M4bg(IwM4E1 < graythresh(IwM4E1)) = 1;

%labeling & dilation
tube1 = false(siz);
L = bwlabeln(M1bg);
tube1(L==1) = true;
SE = strel('disk',16);
mask1 = imdilate(tube1,SE);
mask1 = not(or(mask1,M1bg));

tube2 = false(siz);
L = bwlabeln(M2bg);
tube2(L==1) = true;
mask2 = imdilate(tube2,SE);
mask2 = not(or(mask2,M2bg));

tube3 = false(siz);
L = bwlabeln(M3bg);
tube3(L==1) = true;
mask3 = imdilate(tube3,SE);
mask3 = not(or(mask3,M3bg));

tube4 = false(siz);
L = bwlabeln(M4bg);
tube4(L==1) = true;
mask4 = imdilate(tube4,SE);
mask4 = not(or(mask4,M4bg));
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

imagesc(mask2(:,:,610)');
axis tight equal off
%caxis([0 3])
%colormap(map)
%%
imagesc(pIwM2E2(:,:,610)');
axis tight equal off
caxis([0 0.7])
colormap(gray)
%%
%mask processing
pIwM1E1 = zeros(siz);  pIwM1E1(mask1) = IwM1E1(mask1);   
pIwM1E2 = zeros(siz);  pIwM1E2(mask1) = IwM1E2(mask1);
pIwM1E3 = zeros(siz);  pIwM1E3(mask1) = IwM1E3(mask1);
pIwM1E4 = zeros(siz);  pIwM1E4(mask1) = IwM1E4(mask1);
pIwM1GT = zeros(siz);  pIwM1GT(mask1) = IwM1GT(mask1);
%pIwM1 = zeros(siz);    pIwM1(mask1)   = IwM1(mask1);
mask1 = mask1(:,:,st:en);
pIwM1E1 = pIwM1E1(:,:,st:en); 
pIwM1E2 = pIwM1E2(:,:,st:en); 
pIwM1E3 = pIwM1E3(:,:,st:en); 
pIwM1E4 = pIwM1E4(:,:,st:en);  
%pIwM1 = pIwM1(:,:,st:en);  
pIwM1GT = pIwM1GT(:,:,st:en);

pIwM2E1 = zeros(siz);  pIwM2E1(mask2) = IwM2E1(mask2);   
pIwM2E2 = zeros(siz);  pIwM2E2(mask2) = IwM2E2(mask2);
pIwM2E3 = zeros(siz);  pIwM2E3(mask2) = IwM2E3(mask2);
pIwM2E4 = zeros(siz);  pIwM2E4(mask2) = IwM2E4(mask2);
pIwM2GT = zeros(siz);  pIwM2GT(mask2) = IwM2GT(mask2);
%pIwM2 = zeros(siz);    pIwM2(mask2)   = IwM2(mask2);
mask2 = mask2(:,:,st:en);
pIwM2E1 = pIwM2E1(:,:,st:en); 
pIwM2E2 = pIwM2E2(:,:,st:en); 
pIwM2E3 = pIwM2E3(:,:,st:en); 
pIwM2E4 = pIwM2E4(:,:,st:en); 
%pIwM2 = pIwM2(:,:,st:en);  
pIwM2GT = pIwM2GT(:,:,st:en);

pIwM3E1 = zeros(siz);  pIwM3E1(mask3) = IwM3E1(mask3);   
pIwM3E2 = zeros(siz);  pIwM3E2(mask3) = IwM3E2(mask3);
pIwM3E3 = zeros(siz);  pIwM3E3(mask3) = IwM3E3(mask3);
pIwM3E4 = zeros(siz);  pIwM3E4(mask3) = IwM3E4(mask3);
pIwM3GT = zeros(siz);  pIwM3GT(mask3) = IwM3GT(mask3);
%pIwM3 = zeros(siz);    pIwM3(mask3)   = IwM3(mask3);
mask3 = mask3(:,:,st:en);
pIwM3E1 = pIwM3E1(:,:,st:en); 
pIwM3E2 = pIwM3E2(:,:,st:en); 
pIwM3E3 = pIwM3E3(:,:,st:en); 
pIwM3E4 = pIwM3E4(:,:,st:en);  
%pIwM3 = pIwM3(:,:,st:en);  
pIwM3GT = pIwM3GT(:,:,st:en);

pIwM4E1 = zeros(siz);  pIwM4E1(mask4) = IwM4E1(mask4);   
pIwM4E2 = zeros(siz);  pIwM4E2(mask4) = IwM4E2(mask4);
pIwM4E3 = zeros(siz);  pIwM4E3(mask4) = IwM4E3(mask4);
pIwM4E4 = zeros(siz);  pIwM4E4(mask4) = IwM4E4(mask4);
pIwM4GT = zeros(siz);  pIwM4GT(mask4) = IwM4GT(mask4);
%pIwM4 = zeros(siz);    pIwM4(mask4)   = IwM4(mask4);
mask4 = mask4(:,:,st:en);
pIwM4E1 = pIwM4E1(:,:,st:en); 
pIwM4E2 = pIwM4E2(:,:,st:en); 
pIwM4E3 = pIwM4E3(:,:,st:en); 
pIwM4E4 = pIwM4E4(:,:,st:en);  
%pIwM4 = pIwM4(:,:,st:en);  
pIwM4GT = pIwM4GT(:,:,st:en);

%%
%Save
save_raw(pIwM1E1,'C:\Users\yourb\Desktop\new2\pIwM1E1.raw','*single');
save_raw(pIwM2E1,'C:\Users\yourb\Desktop\new2\pIwM2E1.raw','*single');
save_raw(pIwM3E1,'C:\Users\yourb\Desktop\new2\pIwM3E1.raw','*single');
save_raw(pIwM4E1,'C:\Users\yourb\Desktop\new2\pIwM4E1.raw','*single');

save_raw(pIwM1E2,'C:\Users\yourb\Desktop\new2\pIwM1E2.raw','*single');
save_raw(pIwM2E2,'C:\Users\yourb\Desktop\new2\pIwM2E2.raw','*single');
save_raw(pIwM3E2,'C:\Users\yourb\Desktop\new2\pIwM3E2.raw','*single');
save_raw(pIwM4E2,'C:\Users\yourb\Desktop\new2\pIwM4E2.raw','*single');

save_raw(pIwM1E3,'C:\Users\yourb\Desktop\new2\pIwM1E3.raw','*single');
save_raw(pIwM2E3,'C:\Users\yourb\Desktop\new2\pIwM2E3.raw','*single');
save_raw(pIwM3E3,'C:\Users\yourb\Desktop\new2\pIwM3E3.raw','*single');
save_raw(pIwM4E3,'C:\Users\yourb\Desktop\new2\pIwM4E3.raw','*single');

save_raw(pIwM1E4,'C:\Users\yourb\Desktop\new2\pIwM1E4.raw','*single');
save_raw(pIwM2E4,'C:\Users\yourb\Desktop\new2\pIwM2E4.raw','*single');
save_raw(pIwM3E4,'C:\Users\yourb\Desktop\new2\pIwM3E4.raw','*single');
save_raw(pIwM4E4,'C:\Users\yourb\Desktop\new2\pIwM4E4.raw','*single');

save_raw(pIwM1GT,'C:\Users\yourb\Desktop\new2\pIwM1GT.raw','*uint8');
save_raw(pIwM2GT,'C:\Users\yourb\Desktop\new2\pIwM2GT.raw','*uint8');
save_raw(pIwM3GT,'C:\Users\yourb\Desktop\new2\pIwM3GT.raw','*uint8');
save_raw(pIwM4GT,'C:\Users\yourb\Desktop\new2\pIwM4GT.raw','*uint8');

save_raw(mask1,'C:\Users\yourb\Desktop\new2\pmaskM1.raw','*uint8');
save_raw(mask2,'C:\Users\yourb\Desktop\new2\pmaskM2.raw','*uint8');
save_raw(mask3,'C:\Users\yourb\Desktop\new2\pmaskM3.raw','*uint8');
save_raw(mask4,'C:\Users\yourb\Desktop\new2\pmaskM4.raw','*uint8');
