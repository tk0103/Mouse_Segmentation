%In = pM1E2; InGT = pM1GT; diff1 = 0.25; diff2 = 0.1;
%In = pM2E2; InGT = pM2GT; diff1 = 0.06; diff2 = 0.1;
In = pM3E2; InGT = pM3GT; diff1 = 0.06; diff2 = 0.07;

[N,edges] = histcounts(In(InGT ==1),'BinWidth',0.001);
[~,I] = max(N); modeval(1) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==2),'BinWidth',0.001);
[~,I] = max(N); modeval(2) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==3),'BinWidth',0.001);
[~,I] = max(N); modeval(3) = (edges(I) + edges(I+1)) /2;

[N,edges] = histcounts(In(InGT ==4),'BinWidth',0.001);
[~,I] = max(N); modeval(4) = (edges(I) + edges(I+1)) /2;

Intemp = zeros(siz); mask_blaE2 = zeros(siz); Intemp(InGT == 1) = In(InGT == 1);
mask_blaE2(Intemp > (modeval(1) - diff1) & Intemp < (modeval(1) + diff1) ) = 1; 

Intemp = zeros(siz); mask_LkidE2 = zeros(siz); Intemp(InGT == 2) = In(InGT == 2);
mask_LkidE2(Intemp > (modeval(2) - diff2) & Intemp < (modeval(2) + diff2) ) = 1; 

Intemp = zeros(siz); mask_RkidE2 = zeros(siz); Intemp(InGT == 3) = In(InGT == 3);
mask_RkidE2(Intemp > (modeval(3) - diff2) & Intemp < (modeval(3) + diff2) ) = 1; 

cutM3GT = zeros(siz);  cutM3GT(InGT == 4) = 4;
mask_blaE2 = logical(mask_blaE2); cutM3GT(mask_blaE2) = 1; 
mask_LkidE2 = logical(mask_LkidE2); cutM3GT(mask_LkidE2) = 2;
mask_RkidE2 = logical(mask_RkidE2);cutM3GT(mask_RkidE2) = 3; 
cutM3GT(and(InGT == 1,not(mask_blaE2))) = 5;
cutM3GT(and(InGT == 2,not(mask_LkidE2))) = 6;
cutM3GT(and(InGT == 3,not(mask_RkidE2))) = 7;

cutM3GT = uint8(cutM3GT);
%%
maskM1 = zeros(siz);
maskM1(cutM3GT == 1) = 1;
maskM1(cutM3GT == 2) = 1;
maskM1(cutM3GT == 3) = 1;
%maskM3(cutGTM3 == 4) = 1;
maskM1 = logical(maskM1);
%%
imagesc(cutM3GT(:,:,75)');
    axis tight equal off
    caxis([0 7]);
   % colormap(map)
%%
In = pM1E2; InGT = cutM3GT; Eval = 1; %E2 = 1 E3 = 2 E4 = 3
mu = S.mu; sigma = sqrt(S.Sigma);
%mu = GMMMu; sigma = sqrt(GMMSigma);
edge =[0 0.01:0.01:2.99 3.0]; xlim([0 3.0])
%edge =[0 0.01:0.005:1.295 1.3]; xlim([0 1.3])

hold on
histogram(In(InGT ==1),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==2),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==3),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==4),edge,'Normalization','pdf','EdgeAlpha',0.4);
histogram(In(InGT ==5),edge,'Normalization','pdf','EdgeAlpha',0.4);

y1 = pdf('Normal',edge,mu(1,Eval),sigma(Eval,Eval,1));
plot(edge,y1,'Color',[51 102 255]/255,'LineWidth',2)

y2 = pdf('Normal',edge,mu(2,Eval),sigma(Eval,Eval,2));
plot(edge,y2,'Color',[255 135 0]/255,'LineWidth',2)

y3 = pdf('Normal',edge,mu(3,Eval),sigma(Eval,Eval,3));
plot(edge,y3,'Color',[255 255 0]/255,'LineWidth',2)

y4 = pdf('Normal',edge,mu(4,Eval),sigma(Eval,Eval,4));
plot(edge,y4,'Color',[142 0 204]/255,'LineWidth',2)

y5 = pdf('Normal',edge,mu(5,Eval),sigma(Eval,Eval,5));
plot(edge,y5,'Color',[0 128 0]/255,'LineWidth',2)
%%
tmp1 = pM3E2; tmp2 = pM3E3; tmp3 = pM3E4; mask = cutM3GT;
for k = 1:4
    S.mu(k,1) = mean(tmp1(mask == k));
    S.mu(k,2) = mean(tmp2(mask == k));
    S.mu(k,3) = mean(tmp3(mask == k));
    S.Sigma(:,:,k) = cov(([tmp1(mask == k),tmp2(mask == k),tmp3(mask == k)]));
end
clearvars tmp1 tmp2 tmp3
%%
In = pM3E2; InGT = cutM3GT(newmaskM3);
proval = numel(In(InGT == 5)) / (numel(In(InGT == 2) + numel(In(InGT == 3))));
%%
atlas(:,2) = atlas(:,2) - atlas(:,2) .* proval;
atlas(:,3) = atlas(:,3) - atlas(:,3) .* proval;
atlas(:,5) = atlas(:,2) .* proval + atlas(:,3) .* proval;
%%
[Imap,L,PP,GMMMu,GMMSigma,GMMpro,Feat,likelihood] = AtlasGuidedEM_kubo(Xte,atlas,S,K+1,pmask3,siz,30);
%%
mask = pmask3;
Imap2 = Imap; PPtemp1 = zeros(siz); PPtemp2 = zeros(siz); PPtemp3 = zeros(siz);
PPtemp1(mask) = PP(:,2); PPtemp2(mask) = PP(:,3); PPtemp3(mask) = PP(:,5);

boxval =5;
for i = 1:1:siz(1)
    for j = 1:1:siz(2)
        for k = 1:1:siz(3)
            if(Imap(i,j,k) == 5)
                val1 = sum(PPtemp1(i-boxval:i+boxval,j-boxval:j+boxval,k-boxval:k+boxval));
                val2 = sum(PPtemp2(i-boxval:i+boxval,j-boxval:j+boxval,k-boxval:k+boxval));
                if(val1 > val2)
                    Imap2(i,j,k) = 2;
                    PPtemp1(i,j,k) = PPtemp3(i,j,k);
                    
                else
                    Imap2(i,j,k) = 3;
                    PPtemp2(i,j,k) = PPtemp3(i,j,k);
                end
            end     
        end
    end
end

%%
imagesc(pM1E1(:,:,55)');
axis tight equal off 
%%
JI= CalcuJI(Imap2,pM3GT,K-1);
disp("EM_MAP result")
disp(JI);