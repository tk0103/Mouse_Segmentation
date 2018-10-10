%mask
st = 164; en = 439;
num_of_patch_bla = 4000;
num_of_patch_Lkid = 1000;
num_of_patch_Rkid = 1000;
num_of_patch_bg = 2000;

Label = M1GT;
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
%%
bla(:,4) = 1;
Lkid(:,4) = 5;
Rkid(:,4) = 5;
bg(:,4) = 5;

%%
%Random sampling
idx = randperm(size(bla,1),num_of_patch_bla); bla = bla(idx,:);
idx = randperm(size(Lkid,1),num_of_patch_Lkid); Lkid = Lkid(idx,:);
idx = randperm(size(Rkid,1),num_of_patch_Rkid); Rkid = Rkid(idx,:);
idx = randperm(size(bg,1),num_of_patch_bg); bg = bg(idx,:);
coordi = [bla; Lkid; Rkid; bg;];
%%
%dave
csvwrite("C:\Users\yourb\Desktop\new_tr1.csv",coordi)
%%
imagesc(M3GT(:,:,360)')
axis tight equal
caxis([0 4])
%%
%Test
st = 164; en = 439;
%temp = load_raw('C:\\Users\\yourb\\Desktop\\U_net_result\\7ch\\Results_trM1_ValiM2\\TestResultM3.raw','*uint8');
temp = load_raw('C:\\Users\\yourb\\Desktop\\Result_vali.raw','*uint8');

temp = reshape(temp,siz);
Testresult = zeros(siz); pmask3 = zeros(siz);
pmask3(:,:,st:en) = mask1(:,:,st:en); 
pmask3 = logical(pmask3);
Testresult(pmask3) = temp(pmask3);

%GT
M3GT = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\M2GT.raw','*uint8');
M3GT = reshape(M3GT,siz);

%%
JI = CalcuJI(Testresult,M3GT,4);
disp(JI);

Dice = CalcuDice(Testresult,M3GT,4);
disp(Dice);

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
%%
slice1 = 375;
slice2 = 230;

subplot(2,2,1)
imagesc(Testresult(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,2)
imagesc(M2GT(:,:,slice1)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,3)
imagesc(Testresult(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

subplot(2,2,4)
imagesc(M2GT(:,:,slice2)');
axis tight equal off
caxis([0 4])
colormap(map)

%%
subplot(2,2,1)
imagesc(M1E1(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

subplot(2,2,2)
imagesc(M1E1(:,:,slice1)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

subplot(2,2,3)
imagesc(M1E1(:,:,slice2)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

subplot(2,2,4)
imagesc(M1E1(:,:,slice2)');
axis tight equal off
caxis([0 0.7])
colormap(gray)



%%
temp1 = load_raw('C:\Users\yourb\Documents\GitHub\3D-Unet\Result_test\softou1.raw','*double');
temp2 = load_raw('C:\Users\yourb\Documents\GitHub\3D-Unet\Result_test\softou2.raw','*double');
temp3 = load_raw('C:\Users\yourb\Documents\GitHub\3D-Unet\Result_test\softou3.raw','*double');
temp4 = load_raw('C:\Users\yourb\Documents\GitHub\3D-Unet\Result_test\softou4.raw','*double');
temp1 = reshape(temp1,siz); temp2 = reshape(temp2,siz); temp3 = reshape(temp3,siz); temp4 = reshape(temp4,siz);

%%
slice = 245;
subplot(2,2,1)
imagesc(temp1(:,:,slice)');
axis tight equal off
caxis([0.1 0.3])

subplot(2,2,2);
imagesc(temp2(:,:,slice)');
axis tight equal off
caxis([0.1 0.3])

subplot(2,2,3);
imagesc(temp3(:,:,slice)');
axis tight equal off
caxis([0.1 0.3])

subplot(2,2,4);
imagesc(temp4(:,:,slice)');
axis tight equal off
caxis([0.1 0.3])


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

%%
outout  = (M1E1 - mean(M1E1(M1GT==1))) /  std(M1E1(M1GT==1))*0.8-1;

%%
%histogram
class = 1;
temp1 = outout(M1GT == class);
edges = [-1 -1:0.01:3 3];
hold on
histogram(temp1,edges,'Normalization','probability');
%legend('Mouse1','Mouse2','Mouse3')
