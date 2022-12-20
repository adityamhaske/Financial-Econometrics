function [b_ols,b_bc,b_jk2]=pred_bias_correction(r,x)

n=length(r);
res1=ols(r(2:n),[ones(n-1,1) x(1:(n-1))]);
u_hat=res1.resid;
b_ols=res1.beta(2);
    
res2=ols(x(2:n),[ones(n-1,1) x(1:(n-1))]);
rho_hat=res2.beta(2);
v_hat=res2.resid;
    
phi=u_hat'*v_hat/(v_hat'*v_hat);
    
%stamburgh estimator
b_bc=b_ols+phi*(1+3*rho_hat)/(n-1);
    
%jackknifing (m=2)
    
n1=fix(n/2);
res_jk1=ols(r(2:n1),[ones(n1-1,1) x(1:(n1-1))]);
res_jk2=ols(r((n1+2):n),[ones(n-n1-1,1) x((n1+1):(n-1))]);
    
b_jk2=2*b_ols-(res_jk1.beta(2)+res_jk2.beta(2))/2;