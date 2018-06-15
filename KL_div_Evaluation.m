%%
%KLdiv 1-dimension
ori = pM3E4;
M1B = ori(pM3GT==1);
M1K = ori(or(pM3GT==2,pM3GT==3));
M1E = ori(pM3GT==4);
M1min = min(ori(:));
M1max = max(ori(:));

h = (M1max - M1min)/ (sqrt(size(ori(pmask3),1))/5);
edges = M1min:h:M1max;
[NM1B,~] = histcn(M1B,edges);
[NM1K,~] = histcn(M1K,edges);
[NM1E,~] = histcn(M1E,edges);
NM1B = NM1B ./ sum(NM1B) + eps;
NM1K = NM1K ./ sum(NM1K) + eps;
NM1E = NM1E ./ sum(NM1E) + eps;

P = NM1B;
Q = NM1K;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1B;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1K;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);

%%
%KLdiv 2-dimension
ori1 = pM3E3; ori2 = pM3E4;
M1B = [ori1(pM3GT==1),ori2(pM3GT==1)];
M1K = [ori1(or(pM3GT==2,pM3GT==3)),ori2(or(pM3GT==2,pM3GT==3))];
M1E = [ori1(pM3GT==4),ori2(pM3GT==4)];
M1min1 = min(ori1(:)); M1max1 = max(ori1(:));
M1min2 = min(ori2(:)); M1max2 = max(ori2(:));
h1 = (M1max1 - M1min1)/ (sqrt(size(ori(pmask3),1))/5);
h2 = (M1max2 - M1min2)/ (sqrt(size(ori(pmask3),1))/5);
edges1 = M1min1:h1:M1max1;
edges2 = M1min2:h2:M1max2;
[NM1B,~] = histcn(M1B,edges1,edges2);
[NM1K,~] = histcn(M1K,edges1,edges2);
[NM1E,~] = histcn(M1E,edges1,edges2);
NM1B = NM1B ./ sum(NM1B(:)) + eps;
NM1K = NM1K ./ sum(NM1K(:)) + eps;
NM1E = NM1E ./ sum(NM1E(:)) + eps;

P = NM1B;
Q = NM1K;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1B;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1K;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);

%%
%KLdiv 3-dimension
ori1 = pM3E2; ori2 = pM3E3; ori3 = pM3E4;
label = pM3GT;
M1B = [ori1(label==1),ori2(label==1),ori3(label==1)];
M1K = [ori1(or(label==2,label==3)),ori2(or(label==2,label==3)),ori3(or(label==2,label==3))];
M1E = [ori1(label==4),ori2(label==4),ori3(label==4)];
M1min1 = min(ori1(:)); M1max1 = max(ori1(:));
M1min2 = min(ori2(:)); M1max2 = max(ori2(:));
M1min3 = min(ori3(:)); M1max3 = max(ori3(:));
h1 = (M1max1 - M1min1)/ (sqrt(size(ori(pmask3),1))/5);
h2 = (M1max2 - M1min2)/ (sqrt(size(ori(pmask3),1))/5);
h3 = (M1max3 - M1min3)/ (sqrt(size(ori(pmask3),1))/5);
edges1 = M1min1:h1:M1max1;
edges2 = M1min2:h2:M1max2;
edges3 = M1min3:h3:M1max3;
[NM1B,~] = histcn(M1B,edges1,edges2,edges3);
[NM1K,~] = histcn(M1K,edges1,edges2,edges3);
[NM1E,~] = histcn(M1E,edges1,edges2,edges3);
NM1B = NM1B ./ sum(NM1B(:)) + eps;
NM1K = NM1K ./ sum(NM1K(:)) + eps;
NM1E = NM1E ./ sum(NM1E(:)) + eps;

P = NM1B;
Q = NM1K;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1B;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1K;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);


%%
%KLdiv 4-dimension
ori1 = pM3E1; ori2 = pM3E2; ori3 = pM3E3; ori4 = pM3E4;
label = pM3GT;
M1B = [ori1(label==1),ori2(label==1),ori3(label==1),ori4(label==1)];
M1K = [ori1(or(label==2,label==3)),ori2(or(label==2,label==3)),ori3(or(label==2,label==3)),ori4(or(label==2,label==3))];
M1E = [ori1(label==4),ori2(label==4),ori3(label==4),ori4(label==4)];
M1min1 = min(ori1(:)); M1max1 = max(ori1(:));
M1min2 = min(ori2(:)); M1max2 = max(ori2(:));
M1min3 = min(ori3(:)); M1max3 = max(ori3(:));
M1min4 = min(ori4(:)); M1max4 = max(ori4(:));
h1 = (M1max1 - M1min1)/ (sqrt(size(ori(pmask3),1))/30);
h2 = (M1max2 - M1min2)/ (sqrt(size(ori(pmask3),1))/30);
h3 = (M1max3 - M1min3)/ (sqrt(size(ori(pmask3),1))/30);
h4 = (M1max4 - M1min4)/ (sqrt(size(ori(pmask3),1))/30);
edges1 = M1min1:h1:M1max1;
edges2 = M1min2:h2:M1max2;
edges3 = M1min3:h3:M1max3;
edges4 = M1min4:h4:M1max4;

[NM1B,~] = histcn(M1B,edges1,edges2,edges3,edges4);
[NM1K,~] = histcn(M1K,edges1,edges2,edges3,edges4);
[NM1E,~] = histcn(M1E,edges1,edges2,edges3,edges4);
NM1B = NM1B ./ sum(NM1B(:)) + eps;
NM1K = NM1K ./ sum(NM1K(:)) + eps;
NM1E = NM1E ./ sum(NM1E(:)) + eps;

P = NM1K;
Q = NM1B;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(1,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1B;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(2,1) = (sum(temp1(:)) + sum(temp2(:)))/2;

P = NM1K;
Q = NM1E;
temp1 =  P.*log(P./Q); temp1(isnan(temp1))=0;
temp2 =  Q.*log(Q./P); temp2(isnan(temp2))=0;
dist(3,1) = (sum(temp1(:)) + sum(temp2(:)))/2;
disp(dist);
