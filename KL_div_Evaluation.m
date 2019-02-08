%%
%KLdiv 1-dimension
ori = wM4E1;
GT = M4GT;
mask = mask4;

M1Bla = ori(GT==1);
M1Kid = ori(GT==2);
%M1K = ori(or(GT==2,GT==3));
M1Backg = ori(GT==4);
M1min = min(ori(:));
M1max = max(ori(:));

h = (M1max - M1min)/ (sqrt(size(ori(mask),1))/5);
edges = M1min:h:M1max;
[NM1Bla,~] = histcn(M1Bla,edges);
[NM1Kid,~] = histcn(M1Kid,edges);
[NM1Backg,~] = histcn(M1Backg,edges);
NM1Bla = NM1Bla ./ sum(NM1Bla) + eps;
NM1Kid = NM1Kid ./ sum(NM1Kid) + eps;
NM1Backg = NM1Backg ./ sum(NM1Backg) + eps;

P = NM1Bla;
Q = NM1Kid;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Bla;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Kid;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);

%%
%KLdiv 2-dimension
ori1 = M4E3; 
ori2 = M4E4;
GT =   M4GT;
mask = mask4;

M1Bla = [ori1(GT==1),ori2(GT==1)];
M1Kid = [ori1(GT==2),ori2(GT==2)];
%M1Kid = [ori1(or(GT==2,GT==3)),ori2(or(GT==2,GT==3))];
M1Backg = [ori1(GT==4),ori2(GT==4)];
M1min1 = min(ori1(:)); M1max1 = max(ori1(:));
M1min2 = min(ori2(:)); M1max2 = max(ori2(:));
h1 = (M1max1 - M1min1)/ (sqrt(size(ori(mask),1))/5);
h2 = (M1max2 - M1min2)/ (sqrt(size(ori(mask),1))/5);
edges1 = M1min1:h1:M1max1;
edges2 = M1min2:h2:M1max2;
[NM1Bla,~] = histcn(M1Bla,edges1,edges2);
[NM1Kid,~] = histcn(M1Kid,edges1,edges2);
[NM1Backg,~] = histcn(M1Backg,edges1,edges2);
NM1Bla = NM1Bla ./ sum(NM1Bla(:)) + eps;
NM1Kid = NM1Kid ./ sum(NM1Kid(:)) + eps;
NM1Backg = NM1Backg ./ sum(NM1Backg(:)) + eps;

P = NM1Bla;
Q = NM1Kid;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Bla;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Kid;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);

%%
%KLdiv 3-dimension
ori1 = M4E2; ori2 = M4E3; ori3 = M4E4;
label = M4GT;
mask = mask4;

M1Bla = [ori1(label==1),ori2(label==1),ori3(label==1)];
%M1Kid = [ori1(or(label==2,label==3)),ori2(or(label==2,label==3)),ori3(or(label==2,label==3))];
M1Kid = [ori1(label==2),ori2(label==2),ori3(label==2)];

M1Backg = [ori1(label==4),ori2(label==4),ori3(label==4)];
M1min1 = min(ori1(:)); M1max1 = max(ori1(:));
M1min2 = min(ori2(:)); M1max2 = max(ori2(:));
M1min3 = min(ori3(:)); M1max3 = max(ori3(:));
h1 = (M1max1 - M1min1)/ (sqrt(size(ori(mask),1))/5);
h2 = (M1max2 - M1min2)/ (sqrt(size(ori(mask),1))/5);
h3 = (M1max3 - M1min3)/ (sqrt(size(ori(mask),1))/5);
edges1 = M1min1:h1:M1max1;
edges2 = M1min2:h2:M1max2;
edges3 = M1min3:h3:M1max3;
[NM1Bla,~] = histcn(M1Bla,edges1,edges2,edges3);
[NM1Kid,~] = histcn(M1Kid,edges1,edges2,edges3);
[NM1Backg,~] = histcn(M1Backg,edges1,edges2,edges3);
NM1Bla = NM1Bla ./ sum(NM1Bla(:)) + eps;
NM1Kid = NM1Kid ./ sum(NM1Kid(:)) + eps;
NM1Backg = NM1Backg ./ sum(NM1Backg(:)) + eps;

P = NM1Bla;
Q = NM1Kid;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Bla;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Kid;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);


%%
%KLdiv 4-dimension
ori1 = M4E1; ori2 = M4E2; ori3 = M4E3; ori4 = M4E4;
label = M4GT;
mask = mask4;

M1Bla = [ori1(label==1),ori2(label==1),ori3(label==1),ori4(label==1)];
M1Kid = [ori1(or(label==2,label==3)),ori2(or(label==2,label==3)),ori3(or(label==2,label==3)),ori4(or(label==2,label==3))];
M1Backg = [ori1(label==4),ori2(label==4),ori3(label==4),ori4(label==4)];
M1min1 = min(ori1(:)); M1max1 = max(ori1(:));
M1min2 = min(ori2(:)); M1max2 = max(ori2(:));
M1min3 = min(ori3(:)); M1max3 = max(ori3(:));
M1min4 = min(ori4(:)); M1max4 = max(ori4(:));
h1 = (M1max1 - M1min1)/ (sqrt(size(ori(mask),1))/30);
h2 = (M1max2 - M1min2)/ (sqrt(size(ori(mask),1))/30);
h3 = (M1max3 - M1min3)/ (sqrt(size(ori(mask),1))/30);
h4 = (M1max4 - M1min4)/ (sqrt(size(ori(mask),1))/30);
edges1 = M1min1:h1:M1max1;
edges2 = M1min2:h2:M1max2;
edges3 = M1min3:h3:M1max3;
edges4 = M1min4:h4:M1max4;

[NM1Bla,~] = histcn(M1Bla,edges1,edges2,edges3,edges4);
[NM1Kid,~] = histcn(M1Kid,edges1,edges2,edges3,edges4);
[NM1Backg,~] = histcn(M1Backg,edges1,edges2,edges3,edges4);
NM1Bla = NM1Bla ./ sum(NM1Bla(:)) + eps;
NM1Kid = NM1Kid ./ sum(NM1Kid(:)) + eps;
NM1Backg = NM1Backg ./ sum(NM1Backg(:)) + eps;

P = NM1Kid;
Q = NM1Bla;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Bla;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1Kid;
Q = NM1Backg;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);
