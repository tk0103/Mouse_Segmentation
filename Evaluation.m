%%
%GraphCut Enrgy of each term
seikai = pM1GT(pmask1);
seikai = double(seikai);

Label = CurLabel;
term1 = RP((1:N)'+(Label(:)-1)*N);

Img = pM1E1(pmask1);
Sig = sub2ind(size(Sigmat), Label(GraphModel.Hi), Label(GraphModel.Hj));
Z = Img(GraphModel.Hi)-Img(GraphModel.Hj);
Zval = Z./Sigmat(Sig);
term2 = func(Zval,GraphModel.dist);
term2 = func(-Zval,GraphModel.dist);

shapeval = shape(Sig);
shapeval(shapeval == 2) = E1(shapeval == 2);
shapeval(shapeval == 3) = E2(shapeval == 3);
term3 = shapeval;

term1da = term1(Label == 3);
term2da = term2(Label == 3);
term3da = term3(Label == 3);

disp(sum(term1da(:))*lambda); 
disp(sum(term2da(:))./(1+c));
disp((sum(term3da(:))*c)./(1+c));
term = sum(term1da(:))+sum(term2da(:))+sum(term3da(:));
disp(term);

%%
Label = CurLabel;
E1test = pM1E1(pmask1);
P = GraphModel.Hi(Label(GraphModel.Hi)== 2 & Label(GraphModel.Hj)==4);
Q = GraphModel.Hj(Label(GraphModel.Hi)== 2 & Label(GraphModel.Hj)==4);

sig = 0.0174;
dif = (E1test(P) - E1test(Q))/( sig);
edges = [-2 -2:0.05:3 3];
xlim([-2 3]);
yyaxis left
histogram(dif,edges);
%histogram(dif);

x = edges';
z = exp(-x.^2 /2);
z(x < 0) = 1;

yyaxis right
xlabel('I_p - I_q / \sigma')
plot(x,z,'LineWidth',3);

%%
edge = [0 0:0.01:1.4 1.4];
hold on
histogram(pM1E1(pM1GT ==2),edge,'Normalization','probability','EdgeAlpha',0.4);
histogram(pM1E1(pM1GT ==4),edge,'Normalization','probability','EdgeAlpha',0.4);
xlabel('ŠK’²’l')
ylabel('‘Š‘Î•p“x')
mutest =0.273128082034521;
sigtest = sqrt(0.0014);
y = pdf('Normal',edge,mutest,sigtest)*0.07;
y = y./sum(y(:));
plot(edge,y,'Color',[0.2 0.6 0.2],'LineWidth',2)

%%
mutest =0.273128082034521;
sigtest = sqrt(0.0014);
y = pdf('Normal',edge,mutest,sigtest)*0.07;
y = y./sum(y(:));
plot(edge,y);

%%
%Distance feature
xim = zeros(siz2); yim = zeros(siz2);  zim = zeros(siz2); 
for i = 1:siz2(1)
    xim(i,:,:) = i;
end
for i = 1:siz2(2)
    yim(:,i,:) = i;
end
for i = 1:siz2(3)
    zim(:,:,i) = i;
end

%%
%Calculate kurtosis and skewness
class = 2;
kM1 = [kurtosis(pM1E1(pM1GT == class)); kurtosis(pM1E2(pM1GT == class)); kurtosis(pM1E3(pM1GT == class)); kurtosis(pM1E4(pM1GT == class));];
kM2 = [kurtosis(pM2E1(pM2GT == class)); kurtosis(pM2E2(pM2GT == class)); kurtosis(pM2E3(pM2GT == class)); kurtosis(pM2E4(pM2GT == class));];
kM3 = [kurtosis(pM3E1(pM3GT == class)); kurtosis(pM3E2(pM3GT == class)); kurtosis(pM3E3(pM3GT == class)); kurtosis(pM3E4(pM3GT == class));];

yM1 = [skewness(pM1E1(pM1GT == class)); skewness(pM1E2(pM1GT == class)); skewness(pM1E3(pM1GT == class)); skewness(pM1E4(pM1GT == class));];
yM2 = [skewness(pM2E1(pM2GT == class)); skewness(pM2E2(pM2GT == class)); skewness(pM2E3(pM2GT == class)); skewness(pM2E4(pM2GT == class));];
yM3 = [skewness(pM3E1(pM3GT == class)); skewness(pM3E2(pM3GT == class)); skewness(pM3E3(pM3GT == class)); skewness(pM3E4(pM3GT == class));];


%%
hold on
histogram(pM3E4(pM3GT==1),edge,'Facecolor','y','Normalization','probability','FaceAlpha',0.3);
%histogram(pM3E4(pM3GT==2),edge,'Facecolor','red','Normalization','probability','FaceAlpha',0.3);
%histogram(pM3E4(pM3GT==3),edge,'Facecolor','blue','Normalization','probability','FaceAlpha',0.3);

%%
hold on
scatter(pM3E1(pM3GT==3),pM3E4(pM3GT==3),'.');
scatter(pM3E1(pM3GT==3),pM3E3(pM3GT==3),'.');
scatter(pM3E1(pM3GT==3),pM3E2(pM3GT==3),'.');

scatter(pM3E1(pM3GT==1),pM3E4(pM3GT==1),'.');
scatter(pM3E1(pM3GT==1),pM3E3(pM3GT==1),'.');
scatter(pM3E1(pM3GT==1),pM3E2(pM3GT==1),'.');
axis tight equal
xlim([0.1 2.5]);
ylim([0.1 2.5]);
xlabel('Energy1')
ylabel('Energy2,Energy3,Energy4')
legend('E1&E4(¶t‘Ÿ)','E1&E3(¶t‘Ÿ)','E1&E2(¶t‘Ÿ)','E1&E4(äNã÷)','E1&E3(äNã÷)','E1&E2(äNã÷)','Location','northwest')
%%
class  = 1;
[h1,p1] = lillietest(pM1E3(pM1GT == class));
[h2,p2] = lillietest(pM2E3(pM2GT == class));
[h3,p3] = lillietest(pM3E3(pM3GT == class));

%%
%Energy1
gamma = 1/8;
thval = 0.6;

maxvalE1 = max(M1E1(M1GT ==1));
gammamaskM1 = M1E1>thval; gammamaskM2 = M2E1>thval; gammamaskM3 = M3E1>thval;
gammaoutM1 = M1E1;gammaoutM2 = M2E1;  gammaoutM3 = M3E1;
gammaoutM1(gammamaskM1) = maxvalE1 * (M1E1(gammamaskM1) / maxvalE1).^(1/gamma);
gammaoutM2(gammamaskM2) = maxvalE1 * (M2E1(gammamaskM2) / maxvalE1).^(1/gamma);
gammaoutM3(gammamaskM3) = maxvalE1 * (M3E1(gammamaskM3) / maxvalE1).^(1/gamma);

%%
%Energy2
gamma = 1/8;
thval = 0.7;

maxvalE2 = max(M1E2(M1GT ==1));
gammamaskM1 = M1E2>thval; gammamaskM2 = M2E2>thval; gammamaskM3 = M3E2>thval;
gammaoutM1 = M1E2; gammaoutM2 = M2E2;  gammaoutM3 = M3E2;
gammaoutM1(gammamaskM1) = maxvalE2 * (M1E2(gammamaskM1) / maxvalE2).^(1/gamma);
gammaoutM2(gammamaskM2) = maxvalE2 * (M2E2(gammamaskM2) / maxvalE2).^(1/gamma);
gammaoutM3(gammamaskM3) = maxvalE2 * (M3E2(gammamaskM3) / maxvalE2).^(1/gamma);

%%
%Energy3
gamma = 1/8;
thval = 0.6;

maxvalE3 = max(M1E3(M1GT ==1));
gammamaskM1 = M1E3>thval; gammamaskM2 = M2E3>thval; gammamaskM3 = M3E3>thval;
gammaoutM1 = M1E3; gammaoutM2 = M2E3;  gammaoutM3 = M3E3;
gammaoutM1(gammamaskM1) = maxvalE3 * (M1E3(gammamaskM1) / maxvalE3).^(1/gamma);
gammaoutM2(gammamaskM2) = maxvalE3 * (M2E3(gammamaskM2) / maxvalE3).^(1/gamma);
gammaoutM3(gammamaskM3) = maxvalE3 * (M3E3(gammamaskM3) / maxvalE3).^(1/gamma);

%%
%Energy4
gamma = 1/8;
thval = 0.6;

maxvalE4 = max(M1E4(M1GT ==1));
gammamaskM1 = M1E4>thval; gammamaskM2 = M2E4>thval; gammamaskM3 = M3E4>thval;
gammaoutM1 = M1E4; gammaoutM2 = M2E4;  gammaoutM3 = M3E4;
gammaoutM1(gammamaskM1) = maxvalE4 * (M1E4(gammamaskM1) / maxvalE4).^(1/gamma);
gammaoutM2(gammamaskM2) = maxvalE4 * (M2E4(gammamaskM2) / maxvalE4).^(1/gamma);
gammaoutM3(gammamaskM3) = maxvalE4 * (M3E4(gammamaskM3) / maxvalE4).^(1/gamma);


%%
%histogram
class = 3;

temp1 = gammaoutM1(M1GT == class);
temp2 = gammaoutM2(M2GT == class);
temp3 = gammaoutM3(M3GT == class);

edges = [0 0:0.01:1 1];
hold on
histogram(temp1,edges,'Normalization','probability');
histogram(temp2,edges,'Normalization','probability');
histogram(temp3,edges,'Normalization','probability');
%%
y = prctile(M2E1(M2GT==1),3,1)
%%
edges = [0 0:0.01:3 3];
histogram(M1E1(mask1),edges,'Normalization','probability');

%%
%energyE1 = pM1E1; energyE2 = pM1E2; energyE3 = pM1E3; energyE4 = pM1E4; mask = pmask1;
%energyE1 = pM2E1; energyE2 = pM2E2; energyE3 = pM2E3; energyE4 = pM2E4; mask = pmask2;
energyE1 = pM3E1; energyE2 = pM3E2; energyE3 = pM3E3; energyE4 = pM3E4; mask = pmask3;

IoutM3E1 = (energyE1 - min(energyE1(mask))) ./ (max(energyE1(mask) - min(energyE1(mask)))).*(prctile(energyE1(mask),97,1)- prctile(energyE1(mask),3,1)) + prctile(energyE1(mask),3,1);
IoutM3E2 = (energyE2 - min(energyE2(mask))) ./ (max(energyE2(mask) - min(energyE2(mask)))).*(prctile(energyE2(mask),97,1)- prctile(energyE2(mask),3,1)) + prctile(energyE2(mask),3,1);
IoutM3E3 = (energyE3 - min(energyE3(mask))) ./ (max(energyE3(mask) - min(energyE3(mask)))).*(prctile(energyE3(mask),97,1)- prctile(energyE3(mask),3,1)) + prctile(energyE3(mask),3,1);
IoutM3E4 = (energyE4 - min(energyE4(mask))) ./ (max(energyE4(mask) - min(energyE4(mask)))).*(prctile(energyE4(mask),97,1)- prctile(energyE4(mask),3,1)) + prctile(energyE4(mask),3,1);
%%
%histogram
class = 2;
temp1 = IoutM1E1(pM1GT == class);
temp2 = IoutM2E1(pM2GT == class);
temp3 = IoutM3E1(pM3GT == class);


edges = [0 0:0.01:1.6 1.6];
hold on
histogram(temp1,edges,'Normalization','probability');
histogram(temp2,edges,'Normalization','probability');
histogram(temp3,edges,'Normalization','probability');

%%

%histogram
class = 1;
temp1 = pM1E1(pM1GT == class);
%temp2 = IoutM2E1(pM2GT == class);
temp3 = pM3E1(pM3GT == class);


edges = [0 0:0.005:1.5 1.5];
hold on
histogram(temp1,edges,'Normalization','probability');
%histogram(temp2,edges,'Normalization','probability');
histogram(temp3,edges,'Normalization','probability');

%%
save_raw(IoutM1E1,'C:\\Users\\yourb\\Desktop\\Inputminmax_M3E1.raw','*double')
save_raw(IoutM1E2,'C:\\Users\\yourb\\Desktop\\Inputminmax_M3E2.raw','*double')
save_raw(IoutM1E3,'C:\\Users\\yourb\\Desktop\\Inputminmax_M3E3.raw','*double')
save_raw(IoutM1E4,'C:\\Users\\yourb\\Desktop\\Inputminmax_M3E4.raw','*double')

%%
%Energy1
thval = 0.6;
gammamaskM1 = M1E1>thval; gammamaskM2 = M2E1>thval; gammamaskM3 = M3E1>thval;
gammaoutM1 = M1E1;gammaoutM2 = M2E1;  gammaoutM3 = M3E1;
gammaoutM1(gammamaskM1) = thval;
gammaoutM2(gammamaskM2) = thval;
gammaoutM3(gammamaskM3) = thval;

%%
%Energy4
gamma = 0.1;
thval = prctile(M1E4(M1GT == 2),95,1);
gammamaskM1 = M1E4>=thval; gammamaskM2 = M2E4>=thval; gammamaskM3 = M3E4>=thval;
gammaoutM1E4 = M1E4; gammaoutM2E4 = M2E4;  gammaoutM3E4 = M3E4;
gammaoutM1E4(gammamaskM1) = M1E4(gammamaskM1) .* gamma +(0.9*thval);
gammaoutM2E4(gammamaskM2) = M2E4(gammamaskM2) .* gamma +(0.9*thval);
gammaoutM3E4(gammamaskM3) = M3E4(gammamaskM3) .* gamma +(0.9*thval);
%%
save_raw(IoutM3E1,'C:\\Users\\yourb\\Desktop\\Inputprc99.9_3_M3E1.raw','*double')
save_raw(IoutM3E2,'C:\\Users\\yourb\\Desktop\\Inputprc99.9_3_M3E2.raw','*double')
save_raw(IoutM3E3,'C:\\Users\\yourb\\Desktop\\Inputprc99.9_3_M3E3.raw','*double')
save_raw(IoutM3E4,'C:\\Users\\yourb\\Desktop\\Inputprc99.9_3_M3E4.raw','*double')

%%
edges = [0 0:0.01:2.5 2.5];
hold on
histogram(M1E1(mask1),edges,'Normalization','probability');
histogram(M2E1(mask2),edges,'Normalization','probability');
histogram(M3E1(mask3),edges,'Normalization','probability');

%%
st = 164; en = 439; 
pmask1 = zeros(siz); pmask2 = zeros(siz); pmask3 = zeros(siz);
pmask1(:,:,st:en) = mask1(:,:,st:en); pmask2(:,:,st:en) = mask2(:,:,st:en); pmask3(:,:,st:en) = mask3(:,:,st:en);
pmask1 = logical(pmask1); pmask2 = logical(pmask2); pmask3 = logical(pmask3); 
%%
energyE1 = M1E1; energyE2 = M1E2; energyE3 = M1E3; energyE4 = M1E4; mask = pmask1;
%energyE1 = M2E1; energyE2 = M2E2; energyE3 = M2E3; energyE4 = M2E4; mask = pmask2;
%energyE1 = M3E1; energyE2 = M3E2; energyE3 = M3E3; energyE4 = M3E4; mask = pmask3;

IoutM1E1 = (energyE1 - prctile(energyE1(mask),3,1)) ./ (prctile(energyE1(mask),99.9,1) - prctile(energyE1(mask),3,1));
IoutM1E2 = (energyE2 - prctile(energyE2(mask),3,1)) ./ (prctile(energyE2(mask),99.9,1) - prctile(energyE2(mask),3,1));
IoutM1E3 = (energyE3 - prctile(energyE3(mask),3,1)) ./ (prctile(energyE3(mask),99.9,1) - prctile(energyE3(mask),3,1));
IoutM1E4 = (energyE4 - prctile(energyE4(mask),3,1)) ./ (prctile(energyE4(mask),99.9,1) - prctile(energyE4(mask),3,1));