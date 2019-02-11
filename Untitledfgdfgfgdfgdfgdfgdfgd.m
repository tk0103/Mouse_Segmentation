[sigma,lambda] =ndgrid(0.0021:0.0002:0.004,0.041:0.002:0.07);
lambda = lambda(:);
sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);

%%
[sigma,lambda] =ndgrid(0.0001:0.0004:0.008,0.001:0.002:0.07);
lambda = lambda(:); sigma = sigma(:);
OutputJI = zeros(size(sigma,1),1);