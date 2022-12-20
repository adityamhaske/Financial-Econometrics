function [b_ols,b_boot_bc,t_ols,cv_boot,pv_boot]=boot_naive_2(r,x,level);

%from boot_naive: now add two-sided t-test

%Notes: Only one-predictor. Testing for b=0.

%Input

%r: n by 1, returns
%x: n by 1, predictors
%level: e.g. 0.05

%Output:

%b_ols: 1 by 1
%b_boot_bc: 1 by 1, bias corrected beta (using bootstrap). %bc means bias correction
%t_ols:  1 by 1
%cv_boot: 1 by 3, bootstrap critical values (Left tail, Right tail, both tails)
%pv_boot: 1 by 3, bootstrap p-values (left, right, two-sided)

cv_boot=zeros(1,3);
pv_boot=zeros(1,3);

n=length(r);

bt_m=500; %bootstrap replications

%get the residual to bootstrap
  
res1=ols(r(2:n),[ones(n-1,1) x(1:(n-1))]);   %regression of r
b0=res1.beta(1);
b=res1.beta(2);
b_ols=b;

u0=res1.resid;

%original t-test
t_ols=res1.tstat(2);

res=ols(x(2:end),[ones(n-1,1),x(1:(end-1))]);  %regression of x
rho_0=res.beta(1);
rho=res.beta(2);
ux=res.resid;

u=[u0 ux];

%now bootstrap

b_boot=zeros(1,bt_m);
t_boot=zeros(1,bt_m);

%Now start the bootstrap test.

for bt=1:bt_m;
    %bt
    u_boot=zeros(n-1,2);
    
    for i=1:(n-1); 
        innov_wb=(binornd(1,0.5)-1/2)*2;
        u_boot(i,:)=u(i,:)*innov_wb;        %wild bootstrap
    end;     
    x_boot=zeros(n,1); x_boot(1)=x(1);      
    
    for i=2:n;
        x_boot(i)=rho_0+rho*x_boot(i-1)+u_boot(i-1,2);
    end;
    
    r_boot=zeros(n,1); r_boot(1)=r(1);
    for i=2:n;r_boot(i)=b0+x_boot(i-1,:)*b+u_boot(i-1,1);end;   %null is not imposed
       
    %now perform estimation using the bootstrap data   
    
    res=ols(r_boot(2:end),[ones(n-1,1) x_boot(1:(end-1))]);
    b_boot(bt)=res.beta(2);
    
    %t_boot(bt)=res.tstat(2);  %This is wrong!
    sde=res.beta(2)/res.tstat(2);
    t_boot(bt)=(res.beta(2)-b_ols)/sde;
    
end;

b_boot_bc=b_ols-(mean(b_boot)-b_ols);

cv_boot(1)=quantile(t_boot,level);  %left tail
cv_boot(2)=quantile(t_boot,1-level);  %right tail
cv_boot(3)=quantile(abs(t_boot),1-level);  %two-tailed

su=0;for i=1:bt_m; if (t_boot(i)<t_ols);su=su+1;end;end;pv_boot(1)=su/bt_m;
su=0;for i=1:bt_m; if (t_boot(i)>t_ols);su=su+1;end;end;pv_boot(2)=su/bt_m;
su=0;for i=1:bt_m; if (abs(t_boot(i))>abs(t_ols));su=su+1;end;end;pv_boot(3)=su/bt_m;