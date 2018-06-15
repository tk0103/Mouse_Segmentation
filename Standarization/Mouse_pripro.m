%%
%siz
siz = [544,544,896];
siz2 = [siz(1),siz(2),894];
siz3 = [siz(1),siz(2),864];

%Input_Mouse1
M1E1ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse1\E1.raw','*single');
M1E2ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse1\E2.raw','*single');
M1E3ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse1\E3.raw','*single');
M1E4ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse1\E4.raw','*single');
%M1GTori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse1\GTMouse1.raw','*single');
M1E1ori = reshape(M1E1ori,siz3); M1E2ori = reshape(M1E2ori,siz3); M1E3ori = reshape(M1E3ori,siz3); 
M1E4ori = reshape(M1E4ori,siz3); %M1GTori = reshape(M1GTori,siz3);

%除外
M1E1 = M1E1ori(:,:,3:siz3(3)); M1E2 = M1E2ori(:,:,3:siz3(3)); M1E3 = M1E3ori(:,:,3:siz3(3));
M1E4 = M1E4ori(:,:,3:siz3(3)); %M2GT = M2GTori(:,:,3:siz(3));

%Input_Mouse2
M2E1ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse2\E1.raw','*single');
M2E2ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse2\E2.raw','*single');
M2E3ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse2\E3.raw','*single');
M2E4ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse2\E4.raw','*single');
%M2GTori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse2\GTMouse2.raw','*single');
M2E1ori = reshape(M2E1ori,siz); M2E2ori = reshape(M2E2ori,siz); M2E3ori = reshape(M2E3ori,siz); 
M2E4ori = reshape(M2E4ori,siz); %M2GTori = reshape(M2GTori,siz);

%除外
M2E1 = M2E1ori(:,:,3:siz3(3)); M2E2 = M2E2ori(:,:,3:siz3(3)); M2E3 = M2E3ori(:,:,3:siz3(3));
M2E4 = M2E4ori(:,:,3:siz3(3)); %M2GT = M2GTori(:,:,3:siz(3));

%%
%Input_Mouse3
M3E1ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse3\E1.raw','*single');
M3E2ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse3\E2.raw','*single');
M3E3ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse3\E3.raw','*single');
M3E4ori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse3\E4.raw','*single');
%M3GTori = load_raw('\\tera\user\boku\NZ_data\Mouse\Energy images\Mouse3\GTMouse3.raw','*single');
M3E1ori = reshape(M3E1ori,siz); M3E2ori = reshape(M3E2ori,siz); M3E3ori = reshape(M3E3ori,siz); 
M3E4ori = reshape(M3E4ori,siz); %M3GTori = reshape(M3GTori,siz);

%除外
M3E1 = M3E1ori(:,:,3:siz3(3)); M3E2 = M3E2ori(:,:,3:siz3(3)); M3E3 = M3E3ori(:,:,3:siz3(3));
M3E4 = M3E4ori(:,:,3:siz3(3)); %M3GT = M3GTori(:,:,3:siz(3));

%%
save_raw(M1E1,'C:\Users\yourb\Desktop\M1E1.raw','*single');
save_raw(M1E2,'C:\Users\yourb\Desktop\M1E2.raw','*single');
save_raw(M1E3,'C:\Users\yourb\Desktop\M1E3.raw','*single');
save_raw(M1E4,'C:\Users\yourb\Desktop\M1E4.raw','*single');

save_raw(M2E1,'C:\Users\yourb\Desktop\M2E1.raw','*single');
save_raw(M2E2,'C:\Users\yourb\Desktop\M2E2.raw','*single');
save_raw(M2E3,'C:\Users\yourb\Desktop\M2E3.raw','*single');
save_raw(M2E4,'C:\Users\yourb\Desktop\M2E4.raw','*single');

save_raw(M3E1,'C:\Users\yourb\Desktop\M3E1.raw','*single');
save_raw(M3E2,'C:\Users\yourb\Desktop\M3E2.raw','*single');
save_raw(M3E3,'C:\Users\yourb\Desktop\M3E3.raw','*single');
save_raw(M3E4,'C:\Users\yourb\Desktop\M3E4.raw','*single');

%%
%Input_MD
siz = [544,544,384];
siz3 = [siz(1),siz(2),85:134];

Gd = load_raw('\\tera\user\boku\NZ_data\Multi_Contrast_Phantom\MD images\Gd.raw','*single');
Gold = load_raw('\\tera\user\boku\NZ_data\Multi_Contrast_Phantom\MD images\Gold.raw','*single');
Iodine = load_raw('\\tera\user\boku\NZ_data\Multi_Contrast_Phantom\MD images\Iodine.raw','*single');
HA = load_raw('\\tera\user\boku\NZ_data\Multi_Contrast_Phantom\MD images\HA.raw','*single');
Lipid = load_raw('\\tera\user\boku\NZ_data\Multi_Contrast_Phantom\MD images\Lipid.raw','*single');
Water = load_raw('\\tera\user\boku\NZ_data\Multi_Contrast_Phantom\MD images\Water.raw','*single');

Gd = reshape(Gd,siz); Gold = reshape(Gold,siz); Iodine = reshape(Iodine,siz); 
HA = reshape(HA,siz); Lipid = reshape(Lipid,siz); Water = reshape(Water,siz);
%%
%除外
Gd = Gd(:,:,85:134); Gold = Gold(:,:,85:134); Iodine = Iodine(:,:,85:134);
HA = HA(:,:,85:134); Lipid = Lipid(:,:,85:134); Water = Water(:,:,85:134);

%%
save_raw(Gd,'C:\Users\yourb\Desktop\Gd.raw','*single');
save_raw(Gold,'C:\Users\yourb\Desktop\Gold.raw','*single');
save_raw(Iodine,'C:\Users\yourb\Desktop\Iodine.raw','*single');
save_raw(HA,'C:\Users\yourb\Desktop\HA.raw','*single');
save_raw(Lipid,'C:\Users\yourb\Desktop\Lipid.raw','*single');
save_raw(Water,'C:\Users\yourb\Desktop\Water.raw','*single');