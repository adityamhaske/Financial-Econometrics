function [a_lh,b_lh,R2]=long_horizon_regression_direct(r,x,h);

%from long_horizon_regression_4: only direct estimator.

%version 4: from version 1 (AR(1) model): now add output on intercept, for
%forecasting purpose.

%e.g. [b_lh,b_im,R2]=long_horizon_regression([0;ret_4],dp_4,6)
%x is univariate

%b_lh: direct estimator
%b_im: implied estimator

n=length(r);

%Direct estimator
y=zeros(n-h+1,1);
for i=1:(n-h+1);
    s=0;
    for j=0:(h-1);s=s+r(i+j);end;
    y(i)=s;
end;

res=ols(y(2:(n-h+1)),[ones(n-h,1),x(1:(n-h))]);
a_lh=res.beta(1);
b_lh=res.beta(2);
R2=res.rsqr;