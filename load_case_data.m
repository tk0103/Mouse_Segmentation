%load_path
fileID = fopen('InputPath.txt');
C = textscan(fileID,'%s');
fclose(fileID);
InputPath = C{1,1}; 
InputPath = cell2mat(InputPath);

fileID = fopen('CasePath.txt');
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

siz = [544,544,860];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz);  M1E4 = reshape(M1E4,siz);
M1GT = reshape(M1GT,siz); mask1 = reshape(mask1,siz); mask1 = logical(mask1);
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); M2E4 = reshape(M2E4,siz);
M2GT = reshape(M2GT,siz); mask2 = reshape(mask2,siz); mask2 = logical(mask2);
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); M3E4 = reshape(M3E4,siz); 
M3GT = reshape(M3GT,siz); mask3 = reshape(mask3,siz); mask3 = logical(mask3);
