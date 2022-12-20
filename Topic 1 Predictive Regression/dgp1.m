function [b_ols,b_bc,b_jk2,t_ols]=dgp1(b,rho,del,n,rep);

sig2=[1,del;del,1];

b_ols=zeros(rep,1);
b_bc=zeros(rep,1);
b_jk2=zeros(rep,1);
t_ols=zeros(rep,1);

for kk=1:rep;
    
    r=zeros(n+500,1);x=zeros(n+500,1);
    u=zeros(n+500,2);
    
    for i=1:(n+500);
        u(i,:)=mvnrnd([0,0],sig2);
    end;
    
    x(1)=0;
    for i=2:(n+500);
        x(i)=rho*x(i-1)+u(i,2);
    end;
    
    r(1)=0;
    for i=2:(n+500);
        r(i)=b*x(i-1)+u(i,1);
    end;
    
    x=x(501:n+500,:);
    r=r(501:n+500,:);   %The data are ready
    
    res1=ols(r(2:n),[ones(n-1,1) x(1:(n-1))]);
    u_hat=res1.resid;
    b_ols(kk)=res1.beta(2);
    t_ols(kk)=res1.tstat(2);
    
    res2=ols(x(2:n),[ones(n-1,1) x(1:(n-1))]);
    rho_hat=res2.beta(2);
    v_hat=res2.resid;
    
    phi=u_hat'*v_hat/(v_hat'*v_hat);
    
    b_bc(kk)=b_ols(kk)+phi*(1+3*rho_hat)/(n-1);
    
    %jackknifing (m=2)
    
    n1=fix(n/2);
    res_jk1=ols(r(2:n1),[ones(n1-1,1) x(1:(n1-1))]);
    res_jk2=ols(r((n1+2):n),[ones(n-n1-1,1) x((n1+1):(n-1))]);
    
    b_jk2(kk)=2*b_ols(kk)-(res_jk1.beta(2)+res_jk2.beta(2))/2;
    
end;