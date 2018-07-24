%figure
slice = 200;
subplot(1,4,1);
imagesc(pM3E1(:,:,slice)'); 
axis equal tight off
colormap(gca,'gray');
caxis([0 0.7])

%Colormap
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];

subplot(1,4,2);
imagesc(pM3GT(:,:,slice)'); 
axis equal tight off
colormap(map)
caxis([0 4]);
%caxis([0 0.05])

subplot(1,4,3);

imagesc(Imap(:,:,slice)'); 
axis equal tight off
colormap(map)
caxis([0 4]);
%caxis([0 0.05])

subplot(1,4,4);
imagesc(Output(:,:,slice)'); 
axis equal tight off
colormap(map)
caxis([0 4])

%%
%Voronoi_figure
slice = 270;
imagesc(voronoiOut(:,:,slice));
caxis([0 4])
axis tight equal
%%
map = [0, 0, 0
    parula(4)];

%%
%Save
temp1 = pM1GT;
temp1(temp1 == 4) = 0;

temp2 = Output;
temp2(temp2 == 4) = 0;

save_raw(pM1E1,[InputPath pM1E1 '_Output' '.raw'],'*double');
save_raw(temp1,[InputPath pM1GT '_Output' '.raw'],'*double');
save_raw(temp2,[InputPath Output '_Output' '.raw'],'*double');
