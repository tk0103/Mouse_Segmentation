In = pM3E2; InGT = pM3GT;

modeval(1) = mode(In(InGT == 1));
modeval(2) = mode(In(InGT == 2));
modeval(3) = mode(In(InGT == 3));
modeval(4) = mode(In(InGT == 4));
disp(modeval);

%%
In = pM1E4; InGT = pM1GT;
    histogram(In(InGT == 2),500);
%%
In = pM1E2; InGT = pM1GT; diff1 = 0.25; diff2 = 0.1;

[N,edges] = histcounts(In(InGT ==1),'BinWidth',0.001);
[~,I] = max(N); modeval(1) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==2),'BinWidth',0.001);
[~,I] = max(N); modeval(2) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==3),'BinWidth',0.001);
[~,I] = max(N); modeval(3) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==4),'BinWidth',0.001);
[~,I] = max(N); modeval(4) = (edges(I) + edges(I+1)) /2;

Intemp = zeros(siz2); mask_blaE2 = zeros(siz2); Intemp(InGT == 1) = In(InGT == 1);
mask_blaE2(Intemp > (modeval(1) - diff1) & Intemp < (modeval(1) + diff1) ) = 1; 

Intemp = zeros(siz2); mask_LkidE2 = zeros(siz2); Intemp(InGT == 2) = In(InGT == 2);
mask_LkidE2(Intemp > (modeval(2) - diff2) & Intemp < (modeval(2) + diff2) ) = 1; 

Intemp = zeros(siz2); mask_RkidE2 = zeros(siz2); Intemp(InGT == 3) = In(InGT == 3);
mask_RkidE2(Intemp > (modeval(3) - diff2) & Intemp < (modeval(3) + diff2) ) = 1; 

maskM1 = zeros(siz2); maskM1(InGT == 4) = 4;
maskM1(logical(mask_blaE2)) = 1; maskM1(logical(mask_LkidE2)) = 2; maskM1(logical(mask_RkidE2)) = 3;
%%
for k = 1:4
    tmp1 = pM1E2; tmp2 = pM1E3; tmp3 = pM1E4;
    S.mu(k,1) = mean(tmp1(maskM1 == k));
    S.mu(k,2) = mean(tmp2(maskM1 == k));
    S.mu(k,3) = mean(tmp3(maskM1 == k));
    S.Sigma(:,:,k) = cov(([tmp1(maskM1 == k),tmp2(maskM1 == k),tmp3(maskM1 == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
imagesc(maskM1E3(:,:,160)');
    axis tight equal off
    caxis([0 4]);
    
%%
In = pM1E2; InGT = pM1GT; Eval = 1; %E2 = 1 E3 = 2 E4 = 3
%mu = S.mu; sigma = sqrt(S.Sigma);
mu = GMMMu; sigma = sqrt(GMMSigma);
edge =[0 0.01:0.01:2.99 3.0];
xlim([0 3.0])

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);

y1 = pdf('Normal',edge,mu(1,Eval),sigma(Eval,Eval,1));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

y2 = pdf('Normal',edge,mu(2,Eval),sigma(Eval,Eval,2));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

y3 = pdf('Normal',edge,mu(3,Eval),sigma(Eval,Eval,3));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

y4 = pdf('Normal',edge,mu(4,Eval),sigma(Eval,Eval,4));
plot(edge,y4,'Color',[142 0 204]/255,'LineWidth',2)
%%