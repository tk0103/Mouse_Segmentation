%figure
slice = 191;
subplot(1,4,1);
imagesc(pM1E1(:,:,slice)'); 
axis equal tight off
colormap(gca,'gray');
caxis([0 0.7])

subplot(1,4,2);
imagesc(pM1GT(:,:,slice)'); 
axis equal tight off
colormap(gca,'default')
caxis([0 4]);
%caxis([0 0.05])

subplot(1,4,3);

imagesc(pM2GT(:,:,slice)'); 
axis equal tight off
colormap(gca,'default')
caxis([0 4]);
%caxis([0 0.05])

subplot(1,4,4);
imagesc(pM3GT(:,:,slice)'); 
axis equal tight off
colormap(gca,'default');
caxis([0 4])

%%
%Voronoi_figure
slice = 270;
imagesc(voronoiOut(:,:,slice));
caxis([0 4])
axis tight equal
%%
%Save
temp1 = pM1GT;
temp1(temp1 == 4) = 0;

temp2 = Output;
temp2(temp2 == 4) = 0;

save_raw(pM1E1,[InputPath pM1E1 '_Output' '.raw'],'*double');
save_raw(temp1,[InputPath pM1GT '_Output' '.raw'],'*double');
save_raw(temp2,[InputPath Output '_Output' '.raw'],'*double');
