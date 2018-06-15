%%
slice = 370;
subplot(1,3,1)
imagesc(IwM1E1(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.7])

subplot(1,3,2)
imagesc(IwM2E1(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.7])

subplot(1,3,3)
imagesc(IwM3E1(:,:,slice)');
axis tight equal off
colormap gray
caxis([0 0.7])

%%
slice = 125;
subplot(1,2,1)
imagesc(pIwM1E1(:,:,slice)');
axis tight equal off
colormap default
%caxis([0 3])
colormap gray

subplot(1,2,2)
imagesc(test(:,:,slice)');
axis tight equal off
%caxis([0 3])

slice = 248;
subplot(1,3,3)
imagesc(pIwM2E1(:,:,slice)');
axis tight equal off
%caxis([0 3])
