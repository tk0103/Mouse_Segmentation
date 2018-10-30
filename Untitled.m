temp1 = M1E2 - M1E1;
temp2 = M2E2 - M2E1;
temp3 = M3E2 - M3E1;
%%
hold on
histogram(temp1(mask1));
%%
%%hold on
histogram(temp2(mask2));
%%
histogram(temp3(mask3));
%%
histogram(temp1(M1GT==1));
%%
IoutM3E1 = (energyE1 - prctile(energyE1(mask),3,1)) ./ (prctile(energyE1(mask),99.9,1) - prctile(energyE1(mask),3,1));
IoutM3E2 = (energyE2 - prctile(energyE2(mask),3,1)) ./ (prctile(energyE2(mask),99.9,1) - prctile(energyE2(mask),3,1));
IoutM3E3 = (energyE3 - prctile(energyE3(mask),3,1)) ./ (prctile(energyE3(mask),99.9,1) - prctile(energyE3(mask),3,1));
%%
slice= 230;
subplot(1,2,1);
imagesc(temp1(:,:,slice)');
axis tight equal
%caxis([0 0.1])

subplot(1,2,2);
imagesc(temp2(:,:,slice)');
axis tight equal
%caxis([0 0.1]) 
%%
out = zeros(3,544,544,860);
out(1,:,:,:) = temp3;
out(2,:,:,:) = M3E1;
out(3,:,:,:) = M3E2;
%%
save_raw(out,'C:\\Users\\yourb\\Desktop\\DiffE1E2_M3.raw','*double');
%%
M1GTnew = load_raw('C:\\Users\\yourb\\Desktop\\NZ_unet\\M3GT.raw','*uint8');
M1GTnew = reshape(M1GTnew,[544 544 860]);
%%
slice = 335;
imagesc(M1GTnew(:,:,slice)');
caxis([0 4]);
axis tight equal
%%
M1GTnew(M1GTnew==3) = 2;


%%
save_raw(M1GTnew,'C:\\Users\\yourb\\Desktop\\M3GT3class.raw','*uint8');

%%
temp1 = load_raw('\\tera\user\boku(コスパ改善中)\NZ_data\Mouse\3DVolume\Mouse1\M1E2.raw','*single');
temp2 = load_raw('\\tera\user\boku(コスパ改善中)\NZ_data\Mouse\3DVolume\Mouse2\M2E2.raw','*single');
temp3 = load_raw('\\tera\user\boku(コスパ改善中)\NZ_data\Mouse\3DVolume\Mouse3\M3E2.raw','*single');
%%
temp1 = reshape(temp1,siz);
temp2 = reshape(temp2,siz);
temp3 = reshape(temp3,siz);
%%
imagesc(temp3(:,:,250)');
axis tight equal off
caxis([0 0.7])
colormap(gray)