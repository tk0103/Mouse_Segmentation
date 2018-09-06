%mask
st = 164; en = 439;
num_of_patch_bla = 1000;
num_of_patch_Lkid = 1000;
num_of_patch_Rkid = 1000;
num_of_patch_bg = 2000;

Label = M3GT;
Input_Label = zeros(siz);
Input_Label(:,:,st:en) = Label(:,:,st:en);

[blax,blay,blaz] = ind2sub(siz,find(Input_Label == 1));
[Lkidx,Lkidy,Lkidz] = ind2sub(siz,find(Input_Label == 2));
[Rkidx,Rkidy,Rkidz] = ind2sub(siz,find(Input_Label == 3));
[bgx,bgy,bgz] = ind2sub(siz,find(Input_Label == 4));
bla = [blax,blay,blaz];
Lkid = [Lkidx,Lkidy,Lkidz];
Rkid = [Rkidx,Rkidy,Rkidz];
bg = [bgx,bgy,bgz];

%Random sampling
idx = randperm(size(bla,1),num_of_patch_bla); bla = bla(idx,:);
idx = randperm(size(Lkid,1),num_of_patch_Lkid); Lkid = Lkid(idx,:);
idx = randperm(size(Rkid,1),num_of_patch_Rkid); Rkid = Rkid(idx,:);
idx = randperm(size(bg,1),num_of_patch_bg); bg = bg(idx,:);
coordi = [bla; Lkid; Rkid; bg;];

%dave
csvwrite("C:\Users\yourb\Desktop\tr_M3.csv",coordi)
%%
imagesc(M3GT(:,:,360)')
axis tight equal
caxis([0 4])
%%
%Test
st = 164; en = 439;
temp = load_raw('C:\\Users\\yourb\\Desktop\\U_net_result\\8ch\\Results_trM1_ValiM2_0906\\TestResultM3.raw','*uint8');
temp = reshape(temp,siz);
Testresult = zeros(siz); pmask3 = zeros(siz);
pmask3(:,:,st:en) = mask1(:,:,st:en); pmask3 = logical(pmask3);
Testresult(pmask3) = temp(pmask3);

%GT
M3GT = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\M2GT.raw','*uint8');
M3GT = reshape(M3GT,siz);

%%
%temp = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\8chInput_M3.raw','*double');
temp = load_raw('C:\\Users\\yourb\\Desktop\\8chM3.raw','*double');
%%
temp = reshape(temp,[8 544 544 860]);
%%
temp1 = temp(5,:,:,:); temp1 = squeeze(temp1);
temp2 = temp(6,:,:,:); temp2 = squeeze(temp2);
temp3 = temp(7,:,:,:); temp3 = squeeze(temp3);
temp4 = temp(8,:,:,:); temp4 = squeeze(temp4);
%%
slice = 375;
subplot(2,2,1)
imagesc(temp1(:,:,slice)');
axis tight equal off

subplot(2,2,2);
imagesc(temp2(:,:,slice)');
axis tight equal off
caxis([0 1])

subplot(2,2,3);
imagesc(temp3(:,:,slice)');
axis tight equal off
caxis([0 1])

subplot(2,2,4);
imagesc(temp4(:,:,slice)');
axis tight equal off
caxis([0 1])
%%
JI = CalcuJI(Testresult,M3GT,4);
disp(JI);

Dice = CalcuDice(Testresult,M3GT,4);
disp(Dice);
%%
slice = 407;
imagesc(Testresult(:,:,slice)');
axis tight equal off
%colormap gray
caxis([0 4])
%%
slice = 240;

%Colormap
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

subplot(1,2,1)
imagesc(Testresult(:,:,slice)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(1,2,2)
imagesc(M3GT(:,:,slice)');
axis tight equal off
caxis([0 4])

%%
M3GT = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\M1GT.raw','*uint8');
M3GT = reshape(M3GT,siz);
%%
imagesc(M3GT(:,:,210)');
caxis([0 4])
%%
%histogram
temp1 = ppM1E1(:,:,slice);
temp2 = ppM2E1(:,:,slice);
temp3 = ppM3E1(:,:,slice);
%%
edges = [-2 -2:0.05:1 1];
hold on
histogram(temp1(:),edges);
histogram(temp2(:),edges);
histogram(temp3(:),edges);


%%
%Test_patch creat
stax = 50;
endx = 500;
stay = 50;
endy = 500;
staz = 190;
endz = 439;
out_patch_sizex = 52;
out_patch_sizey = 52;
out_patch_sizez = 52;

num_of_patch = numel(stax:out_patch_sizex:endx) * numel(stay:out_patch_sizey:endy)...
    * numel(staz:out_patch_sizez:endz);
coordi = zeros(num_of_patch,3);
i = 1;

for z = staz:out_patch_sizez:endz
    for y = stay:out_patch_sizey:endy
        for x = stax:out_patch_sizex:endx
            coordi(i,:) = [x y z];
            i = i+1;
        end
    end
end

csvwrite("C:\Users\yourb\Desktop\test_coordinate_nopad52.csv",coordi)









%%
pa_siz = 44/2;
pat_num = 7000;

InPatchCT = zeros(pat_num,4*4*4);
Outresult = zeros(pat_num,4*4*4)+4;
for index = 1:size(coordi,1)
    x = coordi(index,1);
    y = coordi(index,2);
    z = coordi(index,3);
    temp = M1GT(x-pa_siz:x+pa_siz-1,y-pa_siz:y+pa_siz-1,z-pa_siz:z+pa_siz-1);
    temp2 = temp(20:24-1,20:24-1,20:24-1);
    InPatchCT(index,:) = temp2(:);
end

InPatchCT(InPatchCT == 0) = 4;
InPatchCT(InPatchCT == 2) = 1;
InPatchCT(InPatchCT == 3) = 1;

temp1 = InPatchCT == 1;
temp2 = InPatchCT == 4;
class1 = sum(temp1(:));
class2 = sum(temp2(:));

%%
coordinate = csvread("C:\\Users\\yourb\\Documents\\GitHub\\3D-Unet\\Results_trM1_ValiM2\\M2.csv");
GT = zeros(1000,4*4*4);


for index = 1:size(coordinate,1)
    x = coordinate(index,1);
    y = coordinate(index,2);
    z = coordinate(index,3);
    temp = M2GT(x-pa_siz:x+pa_siz-1,y-pa_siz:y+pa_siz-1,z-pa_siz:z+pa_siz-1);
    temp2 = temp(20:24-1,20:24-1,20:24-1);
    GT(index,:) = temp2(:);
end
GT(GT == 0) = 4;
GT(GT == 2) = 1;
GT(GT == 3) = 1;

Result = zeros(1000,4*4*4);
for index = 1:size(coordinate,1)
    x = coordinate(index,1);
    y = coordinate(index,2);
    z = coordinate(index,3);
    temp = Testresult(x-pa_siz:x+pa_siz-1,y-pa_siz:y+pa_siz-1,z-pa_siz:z+pa_siz-1);
    temp2 = temp(20:24-1,20:24-1,20:24-1);
    Result(index,:) = temp2(:);
end
Result(Result == 0) = 4;

%%
testvaliimage = zeros(siz);

temp = Testresult;
temp(temp==0) = 4;
for index = 1:size(coordinate,1)
    x = coordinate(index,1);
    y = coordinate(index,2);
    z = coordinate(index,3);
    temp2 = temp(x-pa_siz:x+pa_siz-1,y-pa_siz:y+pa_siz-1,z-pa_siz:z+pa_siz-1);
    temp2 = temp2(20:24-1,20:24-1,20:24-1);  
    testvaliimage(x-2:x+2-1,y-2:y+2-1,z-2:z+2-1) = temp2;
end
