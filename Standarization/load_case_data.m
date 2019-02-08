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
%%
siz = [544,544,860];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz);  M1E4 = reshape(M1E4,siz);
M1GT = reshape(M1GT,siz);
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); M2E4 = reshape(M2E4,siz);
M2GT = reshape(M2GT,siz);
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); M3E4 = reshape(M3E4,siz); 
M3GT = reshape(M3GT,siz); 
M4E1 = reshape(M4E1,siz); M4E2 = reshape(M4E2,siz); M4E3 = reshape(M4E3,siz); M4E4 = reshape(M4E4,siz); 
M4GT = reshape(M4GT,siz); 

M4E1 = permute(M4E1,[2 1 3]); M4E2 = permute(M4E2,[2 1 3]); M4E3 = permute(M4E3,[2 1 3]);
M4E4 = permute(M4E4,[2 1 3]); M4GT = permute(M4GT,[2 1 3]);
%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

slice = 400;
subplot(1,4,1)
imagesc(M1GT(:,:,242)');
axis tight equal off
colormap(map)
caxis([0 4])

subplot(1,4,2)
imagesc(M2GT(:,:,254)');
axis tight equal off
colormap(map)
caxis([0 4])

subplot(1,4,3)
imagesc(M3GT(:,:,254)');
axis tight equal off
colormap(map)
caxis([0 4])

subplot(1,4,4)
imagesc(M4GT(:,:,248)');
axis tight equal off
colormap(map)
caxis([0 4])


%%
save_raw(M4GT,'C:\Users\yourb\Desktop\label.raw','*uint8');
%%
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
wM1E1 = load_raw([InputPath CasePath{1,:} '.raw'],'*single');
wM2E1 = load_raw([InputPath CasePath{2,:} '.raw'],'*single');
wM3E1 = load_raw([InputPath CasePath{3,:} '.raw'],'*single');
wM4E1 = load_raw([InputPath CasePath{4,:} '.raw'],'*single');

siz = [672 672 792];
wM1E1 = reshape(wM1E1,siz); 
wM2E1 = reshape(wM2E1,siz); 
wM3E1 = reshape(wM3E1,siz); 
siz = [544 544 912];
wM4E1 = reshape(wM4E1,siz);
%%
wM1E1 = wM1E1(65:608,65:608,3:792);
wM2E1 = wM2E1(65:608,65:608,3:792);
wM3E1 = wM3E1(65:608,65:608,3:792);
wM4E1 = wM4E1(:,:,1:790);

siz = [544 544 860];
wM1 = zeros(siz);
wM2 = zeros(siz);
wM3 = zeros(siz);
wM4 = zeros(siz);

wM1(:,:,1:790) = wM1E1;
wM2(:,:,1:790) = wM2E1;
wM3(:,:,1:790) = wM3E1;
wM4(:,:,1:790) = wM4E1;
wM4 = permute(wM4,[2 1 3]);
%%
imagesc(M4E1(:,:,754)');
axis tight equal
