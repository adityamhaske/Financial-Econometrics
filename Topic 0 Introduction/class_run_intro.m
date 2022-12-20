%load 'data1': data on return and dp (n=91)

ret=data1(:,1);
dp=data1(:,2);

subplot(1,2,1)
autocorr(ret)
subplot(1,2,2)
autocorr(dp)

n=length(ret);
res=ols(ret(2:end),[ones(n-1,1) dp(1:end-1)]);
res.beta
res.tstat

[~,~,~,~,stats]=regress(ret(2:end),[ones(n-1,1) dp(1:end-1)])  %The second element is the F-stat. The third element is the p-value of F-stat. F-stat is the square of the t-stat in the one-predictor regression.


res2=hwhite(ret(2:end),[ones(n-1,1) dp(1:end-1)]);
res2.beta    %the same
res2.tstat   %less significant

autocorr(res.resid.*dp(1:end-1))  %no correlation

garchpq(res.resid,1,1)