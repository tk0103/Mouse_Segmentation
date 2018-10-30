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

M1GT = load_raw([InputPath CasePath{13,:} '.raw'],'*uint8');
M2GT = load_raw([InputPath CasePath{14,:} '.raw'],'*uint8');
M3GT = load_raw([InputPath CasePath{15,:} '.raw'],'*uint8');

siz = [544,544,860];
M1E1 = reshape(M1E1,siz); M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz);  M1E4 = reshape(M1E4,siz);
M1GT = reshape(M1GT,siz);
M2E1 = reshape(M2E1,siz); M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); M2E4 = reshape(M2E4,siz);
M2GT = reshape(M2GT,siz);
M3E1 = reshape(M3E1,siz); M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); M3E4 = reshape(M3E4,siz); 
M3GT = reshape(M3GT,siz); 

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

%load_data
M1E1 = load_raw([InputPath CasePath{1,:} '.raw'],'*single');
%M1E2 = load_raw([InputPath CasePath{2,:} '.raw'],'*single');
%M1E3 = load_raw([InputPath CasePath{3,:} '.raw'],'*single');
%M1E4 = load_raw([InputPath CasePath{4,:} '.raw'],'*single');
M2E1 = load_raw([InputPath CasePath{5,:} '.raw'],'*single');
M2E2 = load_raw([InputPath CasePath{6,:} '.raw'],'*single');
%M2E3 = load_raw([InputPath CasePath{7,:} '.raw'],'*single');
%M2E4 = load_raw([InputPath CasePath{8,:} '.raw'],'*single');
M3E1 = load_raw([InputPath CasePath{9,:} '.raw'],'*single');
%M3E2 = load_raw([InputPath CasePath{10,:} '.raw'],'*single');
%M3E3 = load_raw([InputPath CasePath{11,:} '.raw'],'*single');
%M3E4 = load_raw([InputPath CasePath{12,:} '.raw'],'*single');

siz = [672 672 792];
M1E1 = reshape(M1E1,siz); 
%M1E2 = reshape(M1E2,siz); M1E3 = reshape(M1E3,siz);  M1E4 = reshape(M1E4,siz);
M2E1 = reshape(M2E1,siz); 
%M2E2 = reshape(M2E2,siz); M2E3 = reshape(M2E3,siz); M2E4 = reshape(M2E4,siz);
M3E1 = reshape(M3E1,siz); 
%M3E2 = reshape(M3E2,siz); M3E3 = reshape(M3E3,siz); M3E4 = reshape(M3E4,siz); 

wM1E1 = M1E1(65:608,65:608,3:792);
% M1E2 = M1E2(65:608,65:608,3:792);
% M1E3 = M1E3(65:608,65:608,3:792);
% M1E4 = M1E4(65:608,65:608,3:792);
wM2E1 = M2E1(65:608,65:608,3:792);
% M2E2 = M2E2(65:608,65:608,3:792);
% M2E3 = M2E3(65:608,65:608,3:792);
% M2E4 = M2E4(65:608,65:608,3:792);
wM3E1 = M3E1(65:608,65:608,3:792);
% M3E2 = M3E2(65:608,65:608,3:792);
% M3E3 = M3E3(65:608,65:608,3:792);
% M3E4 = M3E4(65:608,65:608,3:792);
