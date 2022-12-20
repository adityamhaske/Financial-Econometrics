function [MSE_0,MSE_hm,MSE,R2,yplot_true,yplot_hm,yplot_direct]=oos_rsq_m524(r,x,h,f_start);

%from oos_rsq_5: Now only the direct estimator (for the class M524). An
%important error is corrected (regarding which predictor value should be used)

%f_start is a fraction 0<f_start<1 (where to start forecast).

%Compute out-of-sample R-square over horizon h

n=length(r);
n1=fix(n*f_start);  %estimation sample
%n1=fix(n/2);
n2=n-n1;

%aggregate dependent variable (containing true values)
y=zeros(n-h+1,1);
for i=1:(n-h+1);
    s=0;
    for j=0:(h-1);s=s+r(i+j);end;
    y(i)=s;
end;

y_true=y((n1+1):(n-h+1));

%historical mean
y_hm=zeros(n2-h+1,1);
for i=1:length(y_hm);
    y_hm(i)=mean(y(1:(n1-h+i)));  %expanding window
end;

MSE_hm=mean((y_hm-y_true).^2);
MSE_0=mean((0-y_true).^2);

%direct method
y_direct=zeros(n2-h+1,1);
for i=1:length(y_direct);
    [a_lh,b_lh,~]=long_horizon_regression_direct(r(1:(n1-h+i)),x(1:(n1-h+i)),h);
    %y_direct(i)=a_lh+b_lh*x(n1-h+i);  %wrong here
    y_direct(i)=a_lh+b_lh*x(n1-1+i);  
end;

MSE=mean((y_direct-y_true).^2);

R2=1-MSE/MSE_hm;

yplot_true=[y(1:n1);y_true];  %These are for plotting purpose
yplot_hm=[y(1:n1);y_hm];
yplot_direct=[y(1:n1);y_direct];