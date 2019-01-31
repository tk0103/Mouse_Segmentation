bn= [0.6551,0.8797,0.7259,	0.8988]';
bw=	[0.7664,0.8538,	0.8682,	0.3976]';
ln=	[0.6288,0.7265,	0.6719,	0.796]';
lw= [0.5277,0.5981,	0.4485,	0.3426]';
rn= [0.7432,0.7785,	0.6519,	0.5236]';
rw= [0.708,	0.6382,	0.2884,	0.3354]';

wn=[0.6757	,0.7949,	0.683233333,	0.739466667];
ww=[0.667366667,	0.6967	,0.535033333,0.358533333];

%%
boxplot([bn,bw,ln,lw,rn,rw],'Labels',{'Bladder(Narrow)','Bladder(Wide)'...
    ,'L.Kidney(Narrow)','L.Kidney(Wide)','R.Kidney(Narrow)','R.Kidney(Wide)'});

%%
[p1,h1] = ranksum(bn,bw)
%%
[p,h] = ranksum(ln,lw)
%%
[p1,h1] = ranksum(rn,rw)

%%
[p1,h1] = ranksum(wn,ww)

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
out = Imap2;
out(Imap2 == 4) = 0;
%imagesc(out(110:290,150:330,200)');
%imagesc(out(260:370,250:360,73)');
imagesc(out(:,:,200)');

axis tight equal off
colormap(map)
caxis([0 4])

%%
imagesc(wM3E1(260:370,250:360,73)');
axis tight equal off
caxis([0 0.7])
colormap(gray)

%%
map = [0, 0, 0
    0.1, 0.5, 0.8
    0.2, 0.7, 0.6
    0.8, 0.7, 0.3
    0.9, 0.9, 0];
out = Imap2;
out(Imap2 == 4) = 0;
imagesc(out(100:290,150:340,200)');
%imagesc(out(:,:,200)');
axis tight equal off
caxis([0 4])
colormap(map);