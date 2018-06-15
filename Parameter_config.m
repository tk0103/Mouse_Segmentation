%clip Input image
st = 164; en = 439;

%the number of class
K=4;

%atlas mouse1
sig1 = 5; %bladder
sig2 = 3; %kidneys

%GraphCut
lambda = 0.05;
h = 0.5;
c=0;

%GraphCut prior
graydiff = zeros(K);
graydiff(1,2) = 1; graydiff(1,3) = 1; graydiff(2,1) = 1; graydiff(3,1) = 1; graydiff(2,3) = 1;graydiff(3,2) = 1;
graydiff(1,4) = 2; graydiff(2,4) = 2; graydiff(3,4) = 2;
graydiff(4,1) = 3; graydiff(4,2) = 3; graydiff(4,3) = 3;

shape = zeros(K);
shape(1,2) = 1; shape(1,3) = 1; shape(2,1) = 1; shape(3,1) = 1; shape(2,3) = 1;shape(3,2) = 1;
shape(1,4) = 0; shape(4,1) = 0;
shape(2,4) = 2; shape(4,2) = 2;
shape(3,4) = 3; shape(4,3) = 3;


%%
%GraphCut Gridserch
[lambda,h,c] =ndgrid(0.01:0.01:0.04,0.1:0.2:1,0.2:0.2:1.0);
lambda = lambda(:);
h = h(:);
c = c(:);
sumJI = zeros(size(h,1),1);