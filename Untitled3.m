maxcomp = zeros(siz); 
L1 = bwconncomp(Imap2 == 3);
[~,idx] = max(cellfun(@numel,L1.PixelIdxList)); 
maxcomp(L1.PixelIdxList{idx}) = 1;
maxcomp = logical(maxcomp);
[XX,YY,ZZ] = meshgrid(1:siz(1),1:siz(2),1:siz(3));
out = [XX(maxcomp) YY(maxcomp) ZZ(maxcomp)];
[coeff,~,latent] = pca(out,'NumComponents',3);

original = zeros(siz);
r = ( ((XX-272)./ (sqrt(latent(1)))).^2 ...
    + ((YY-272)./ (sqrt(latent(2)))).^2 ...
    + ((ZZ-100)./ (sqrt(latent(3)))).^2 );
original( r <= 1 ) = 1;
original = signed_EUDT(original);

mo = mean(out,1); move = eye(4); move2 = eye(4); rotate = eye(4); 
rotate(1:3,1:3) = coeff;
move2(1,4) = mo(1); move2(2,4) = mo(2); move2(3,4) = mo(3);
move(1,4) = -272; move(2,4) = -272; move(3,4) = -100;
invmat = inv(move2*rotate*move);

Ax = single(1:siz(2));
Ay = single(1:siz(1));
Az = single(1:siz(3));
Ax = reshape(Ax,1,siz(1),1);
Ay = reshape(Ay,siz(1),1,1);
Az = reshape(Az,1,1,siz(3));

Gx = (Ax * invmat(1,1) + Ay * invmat(1,2) + Az * invmat(1,3)) + invmat(1,4);
Gy = (Ax * invmat(2,1) + Ay * invmat(2,2) + Az * invmat(2,3)) + invmat(2,4);
Gz = (Ax * invmat(3,1) + Ay * invmat(3,2) + Az * invmat(3,3)) + invmat(3,4);

Re = zeros(siz);
F = griddedInterpolant(original,'linear','nearest');

Iw = F(Gx,Gy,Gz);
Re(Iw <=0 ) =1;

clearvars Ax Ay Az Gx Gy Gz maxcomp XX YY ZZ invmat
%%
imagesc(Imap2(:,:,200)');
axis tight equal
%%
imagesc(Iw(:,:,200)');
axis tight equal
caxis([-10 50])