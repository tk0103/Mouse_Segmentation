%narrow
%load_path
fileID = fopen('InputPath_narrow.txt');
C = textscan(fileID,'%s');
fclose(fileID);
InputPath = C{1,1}; 
InputPath = cell2mat(InputPath);

fileID = fopen('CasePath_narrow.txt');
C = textscan(fileID,'%s');
fclose(fileID);
CasePath = C{1,1};
%%
%load_data
M1E1 = load_raw([InputPath CasePath{1,:} '.raw'],'*double');
M1E2 = load_raw([InputPath CasePath{2,:} '.raw'],'*double');
M1E3 = load_raw([InputPath CasePath{3,:} '.raw'],'*double');
M1E4 = load_raw([InputPath CasePath{4,:} '.raw'],'*double');
M2E1 = load_raw([InputPath CasePath{5,:} '.raw'],'*double');
M2E2 = load_raw([InputPath CasePath{6,:} '.raw'],'*double');
M2E3 = load_raw([InputPath CasePath{7,:} '.raw'],'*double');
M2E4 = load_raw([InputPath CasePath{8,:} '.raw'],'*double');
M3E1 = load_raw([InputPath CasePath{9,:} '.raw'],'*double');
M3E2 = load_raw([InputPath CasePath{10,:} '.raw'],'*double');
M3E3 = load_raw([InputPath CasePath{11,:} '.raw'],'*double');
M3E4 = load_raw([InputPath CasePath{12,:} '.raw'],'*double');

M1GT = load_raw([InputPath CasePath{13,:} '.raw'],'*uint8');
M2GT = load_raw([InputPath CasePath{14,:} '.raw'],'*uint8');
M3GT = load_raw([InputPath CasePath{15,:} '.raw'],'*uint8');

mask1 = load_raw([InputPath CasePath{16,:} '.raw'],'*uint8');
mask2 = load_raw([InputPath CasePath{17,:} '.raw'],'*uint8');
mask3 = load_raw([InputPath CasePath{18,:} '.raw'],'*uint8');

siz = [544 544  860];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz);  M1E4 = reshape(M1E4,siz);
M1GT = reshape(M1GT,siz); mask1 = reshape(mask1,siz); mask1 = logical(mask1);
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); M2E4 = reshape(M2E4,siz);
M2GT = reshape(M2GT,siz); mask2 = reshape(mask2,siz); mask2 = logical(mask2);
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); M3E4 = reshape(M3E4,siz); 
M3GT = reshape(M3GT,siz); mask3 = reshape(mask3,siz); mask3 = logical(mask3);

%%
imagesc(M1E1(:,:,272)');
colormap(gray);

%%
%Wide
%load_path
fileID = fopen('InputPath_wide.txt');
C = textscan(fileID,'%s');
fclose(fileID);
InputPath = C{1,1}; 
InputPath = cell2mat(InputPath);

fileID = fopen('CasePath_wide.txt');
C = textscan(fileID,'%s');
fclose(fileID);
CasePath = C{1,1};
%%
%load_data
wM1E1 = load_raw([InputPath CasePath{1,:} '.raw'],'*double');
wM2E1 = load_raw([InputPath CasePath{2,:} '.raw'],'*double');
wM3E1 = load_raw([InputPath CasePath{3,:} '.raw'],'*double');

wM1GT = load_raw([InputPath CasePath{4,:} '.raw'],'*uint8');
wM2GT = load_raw([InputPath CasePath{5,:} '.raw'],'*uint8');
wM3GT = load_raw([InputPath CasePath{6,:} '.raw'],'*uint8');

wmask1 = load_raw([InputPath CasePath{7,:} '.raw'],'*uint8');
wmask2 = load_raw([InputPath CasePath{8,:} '.raw'],'*uint8');
wmask3 = load_raw([InputPath CasePath{9,:} '.raw'],'*uint8');

siz2 = [544,544,790];
wM1E1 = reshape(wM1E1,siz2); wM1GT = reshape(wM1GT,siz); wmask1 = reshape(wmask1,siz); wmask1 = logical(wmask1);
wM2E1 = reshape(wM2E1,siz2); wM2GT = reshape(wM2GT,siz); wmask2 = reshape(wmask2,siz); wmask2 = logical(wmask2);
wM3E1 = reshape(wM3E1,siz2); wM3GT = reshape(wM3GT,siz); wmask3 = reshape(wmask3,siz); wmask3 = logical(wmask3);

%%
temp = zeros([4 544 544 860]);
temp(1,:,:,:) = M3GT;
temp(2,:,:,:) = M3GT;
temp(3,:,:,:) = M3GT;
temp(4,:,:,:) = M3GT;
%%
save_raw(temp,'C:\\Users\\yourb\\Desktop\\4chLabel_M31.raw','*uint8');




%%
%normilizatoon & atlas data
PP1 = zeros(siz2); PP2 = zeros(siz2); PP3 = zeros(siz2);  

PP1(pmask2) = PP(:,1);
PP2(pmask2) = PP(:,2);
PP3(pmask2) = PP(:,3);

PP1 = imgaussfilt3(PP1,5);
PP2 = imgaussfilt3(PP2,5);
PP3 = imgaussfilt3(PP3,5);
%%
newPP1 = zeros(siz); newPP2 = zeros(siz); newPP3 = zeros(siz); newPP4 = zeros(siz);
newPP1(:,:,st:en) = PP1;
newPP2(:,:,st:en) = PP2;
newPP3(:,:,st:en) = PP3;
newPP4(:,:,st:en) = PP4;
%%
maxval = max(pM1E1(:)); minval = min(pM1E1(:));
normM1E1 = (M1E1 - minval)/(maxval - minval);
maxval = max(pM1E2(:)); minval = min(pM1E2(:));
normM1E2 = (M1E2 - minval)/(maxval - minval);
maxval = max(pM1E3(:)); minval = min(pM1E3(:));
normM1E3 = (M1E3 - minval)/(maxval - minval);
maxval = max(pM1E4(:)); minval = min(pM1E4(:));
normM1E4 = (M1E4 - minval)/(maxval - minval);

maxval = max(pM2E1(:)); minval = min(pM2E1(:));
normM2E1 = (M2E1 - minval)/(maxval - minval);
maxval = max(pM2E2(:)); minval = min(pM2E2(:));
normM2E2 = (M2E2 - minval)/(maxval - minval);
maxval = max(pM2E3(:)); minval = min(pM2E3(:));
normM2E3 = (M2E3 - minval)/(maxval - minval);
maxval = max(pM2E4(:)); minval = min(pM2E4(:));
normM2E4 = (M2E4 - minval)/(maxval - minval);

maxval = max(pM3E1(:)); minval = min(pM3E1(:));
normM3E1 = (M3E1 - minval)/(maxval - minval);
maxval = max(pM3E2(:)); minval = min(pM3E2(:));
normM3E2 = (M3E2 - minval)/(maxval - minval);
maxval = max(pM3E3(:)); minval = min(pM3E3(:));
normM3E3 = (M3E3 - minval)/(maxval - minval);
maxval = max(pM3E4(:)); minval = min(pM3E4(:));
normM3E4 = (M3E4 - minval)/(maxval - minval);

%%
output = zeros([6 544 544 860]);
output(1,:,:,:) = M2E1;
output(2,:,:,:) = M2E2;
output(3,:,:,:) = M2E3;
output(4,:,:,:) = M2E4;
output(5,:,:,:) = newPP2;
output(6,:,:,:) = newPP3;

save_raw(output,'C:\\Users\\yourb\\Desktop\\6chInput_M2_notnorm.raw','*double')

%%
M18ch = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\7chInput_M1_notnorm.raw','*double');
M28ch = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\7chInput_M2_notnorm.raw','*double');
M38ch = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\7chInput_M3_notnorm.raw','*double');
%%
siz = [7 544 544 860];
M18ch = reshape(M18ch,siz); M28ch = reshape(M28ch,siz); M38ch = reshape(M38ch,siz);


%%
tempM1 = zeros([6 544  544 860]);
tempM2 = zeros([6 544  544 860]);
tempM3 = zeros([6 544  544 860]);
%%
tempM1(1:4,:,:,:) = M18ch(1:4,:,:,:);
tempM2(1:4,:,:,:) = M28ch(1:4,:,:,:);
tempM3(1:4,:,:,:) = M38ch(1:4,:,:,:);
%%
tempM1(5:6,:,:,:) = M18ch(6:7,:,:,:);
tempM2(5:6,:,:,:) = M28ch(6:7,:,:,:);
tempM3(5:6,:,:,:) = M38ch(6:7,:,:,:);
%%
save_raw(tempM1,'C:\\Users\\yourb\\Desktop\\6chInput_M1_notnorm.raw','*double')
save_raw(tempM2,'C:\\Users\\yourb\\Desktop\\6chInput_M2_notnorm.raw','*double')
save_raw(tempM3,'C:\\Users\\yourb\\Desktop\\6chInput_M3_notnorm.raw','*double')

%%
M28ch = load_raw('C:\\Users\\yourb\\Desktop\\6chInput_M2.raw','*double');
M28ch = reshape(M28ch,[6 544 544 860]);
%%
temp1 = M28ch(1,:,:,:);
temp2 = M28ch(5,:,:,:);
temp3 = M28ch(6,:,:,:);



%%
temp1 = squeeze(temp1);
temp2 = squeeze(temp2);
temp3 = squeeze(temp3);
%%
slice = 370;
subplot(2,2,1)
imagesc(temp1(:,:,slice)');
axis tight equal

subplot(2,2,2)
imagesc(temp2(:,:,slice)');
axis tight equal
caxis([0 1])

subplot(2,2,3)
imagesc(temp3(:,:,slice)');
axis tight equal
caxis([0 1])

%%
softout1 = load_raw('C:\\Users\\yourb\\Desktop\\softou1.raw','*double');
softout2 = load_raw('C:\\Users\\yourb\\Desktop\\softou2.raw','*double');
softout3 = load_raw('C:\\Users\\yourb\\Desktop\\softou3.raw','*double');
softout4 = load_raw('C:\\Users\\yourb\\Desktop\\softou4.raw','*double');

%%
siz = [544 544 860];
softout1 = reshape(softout1,siz);
softout2 = reshape(softout2,siz);
softout3 = reshape(softout3,siz);
softout4 = reshape(softout4,siz);

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

slice=  220;
subplot(2,2,1)
imagesc(softout1(:,:,slice)');
axis tight equal off
%caxis([0.4 0.6])
caxis([0.46 0.47])
colormap(gray)

subplot(2,2,2)
imagesc(softout2(:,:,slice)');
axis tight equal off
caxis([0.17 0.18])
%caxis([0 0.5])
colormap(map)

subplot(2,2,3)
imagesc(softout3(:,:,slice)');
axis tight equal off
caxis([0 0.5])

subplot(2,2,4)
imagesc(softout4(:,:,slice)');
axis tight equal off
caxis([0 0.5])
%%
