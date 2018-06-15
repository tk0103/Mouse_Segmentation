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