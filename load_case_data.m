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
M1E1 = load_raw([InputPath CasePath{1,:} '.raw'],'*single');
M1E2 = load_raw([InputPath CasePath{2,:} '.raw'],'*single');
M1E3 = load_raw([InputPath CasePath{3,:} '.raw'],'*single');
M1E4 = load_raw([InputPath CasePath{4,:} '.raw'],'*single');
M2E1 = load_raw([InputPath CasePath{5,:} '.raw'],'*single');
M2E2 = load_raw([InputPath CasePath{6,:} '.raw'],'*single');
M2E3 = load_raw([InputPath CasePath{7,:} '.raw'],'*single');
M2E4 = load_raw([InputPath CasePath{8,:} '.raw'],'*single');
M3E1 = load_raw([InputPath CasePath{9,:} '.raw'],'*single');
M3E2 = load_raw([InputPath CasePath{10,:} '.raw'],'*single');
M3E3 = load_raw([InputPath CasePath{11,:} '.raw'],'*single');
M3E4 = load_raw([InputPath CasePath{12,:} '.raw'],'*single');
M4E1 = load_raw([InputPath CasePath{13,:} '.raw'],'*single');
M4E2 = load_raw([InputPath CasePath{14,:} '.raw'],'*single');
M4E3 = load_raw([InputPath CasePath{15,:} '.raw'],'*single');
M4E4 = load_raw([InputPath CasePath{16,:} '.raw'],'*single');

M1GT = load_raw([InputPath CasePath{17,:} '.raw'],'*uint8');
M2GT = load_raw([InputPath CasePath{18,:} '.raw'],'*uint8');
M3GT = load_raw([InputPath CasePath{19,:} '.raw'],'*uint8');
M4GT = load_raw([InputPath CasePath{20,:} '.raw'],'*uint8');

mask1 = load_raw([InputPath CasePath{21,:} '.raw'],'*uint8');
mask2 = load_raw([InputPath CasePath{22,:} '.raw'],'*uint8');
mask3 = load_raw([InputPath CasePath{23,:} '.raw'],'*uint8');
mask4 = load_raw([InputPath CasePath{24,:} '.raw'],'*uint8');
%%
siz = [544 544  287];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz);  M1E4 = reshape(M1E4,siz);
M1GT = reshape(M1GT,siz); mask1 = reshape(mask1,siz); mask1 = logical(mask1);
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); M2E4 = reshape(M2E4,siz);
M2GT = reshape(M2GT,siz); mask2 = reshape(mask2,siz); mask2 = logical(mask2);
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); M3E4 = reshape(M3E4,siz); 
M3GT = reshape(M3GT,siz); mask3 = reshape(mask3,siz); mask3 = logical(mask3);
M4E1 = reshape(M4E1,siz); M4E2 = reshape(M4E2,siz); M4E3 = reshape(M4E3,siz); M4E4 = reshape(M4E4,siz); 
M4GT = reshape(M4GT,siz); mask4 = reshape(mask4,siz); mask4 = logical(mask4);
%%
imagesc(M1GT(:,:,226)');
%colormap(gray);

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

siz3 = [544,544,790];
wM1E1 = reshape(wM1E1,siz3); 
wM2E1 = reshape(wM2E1,siz3);
wM3E1 = reshape(wM3E1,siz3); 
wM1GT = reshape(wM1GT,siz);
wM2GT = reshape(wM2GT,siz); 
wM3GT = reshape(wM3GT,siz); 
wmask1 = reshape(wmask1,siz); wmask1 = logical(wmask1);
wmask2 = reshape(wmask2,siz); wmask2 = logical(wmask2);
wmask3 = reshape(wmask3,siz); wmask3 = logical(wmask3);
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
output = zeros([7 544 544 860]);
output(1,:,:,:) = IoutM3E1;
output(2,:,:,:) = IoutM3E2;
output(3,:,:,:) = IoutM3E3;
output(4,:,:,:) = IoutM3E4;
output(5,:,:,:) = PP1;
output(6,:,:,:) = PP2;
output(7,:,:,:) = PP3;
%%
save_raw(output,'C:\\Users\\yourb\\Desktop\\7chInput_M3_prc99.93.raw','*double')

%%
M17ch = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\7chInput_M1_prc99.93.raw','*double');
M27ch = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\7chInput_M2_prc99.93.raw','*double');
M37ch = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\7chInput_M3_prc99.93.raw','*double');
siz = [7 544 544 860];
M17ch = reshape(M17ch,siz); M27ch = reshape(M27ch,siz); M37ch = reshape(M37ch,siz);

%%
PP1 = M37ch(5,:,:,:);
PP2 = M37ch(6,:,:,:);
PP3 = M37ch(7,:,:,:);
PP1 = squeeze(PP1);
PP2 = squeeze(PP2);
PP3 = squeeze(PP3);

%%
M3E1 = M37ch(1,:,:,:);
M3E2 = M37ch(2,:,:,:);
M3E3 = M37ch(3,:,:,:);
M3E4 = M37ch(4,:,:,:);
M3E1 = squeeze(M3E1);
M3E2 = squeeze(M3E2);
M3E3 = squeeze(M3E3);
M3E4 = squeeze(M3E4);
%%
slice1 = 230;
slice2  = 365;
subplot(1,3,1)
imagesc(PP1(:,:,slice1)');
axis tight equal off
caxis([0 1])
colormap(map)

subplot(1,3,2)
imagesc(PP2(:,:,slice2)');
axis tight equal off
%caxis([0 1])

subplot(1,3,3)
imagesc(PP3(:,:,slice2)');
axis tight equal off
%caxis([0 1])
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];


%%
out = M3GT;
out(out==4) = 0;
%%
slice=  230;
imagesc(out(:,:,slice)');

axis tight equal off
caxis([0 4])
colormap(map)
%%
imagesc(pwM1E1(:,:,slice2)');
colormap(gray)
caxis([0 0.7])
axis tight equal off
%%
 rectangle('Position',[315,250,70,70],'FaceColor','none','EdgeColor','r',...
    'LineWidth',1)
