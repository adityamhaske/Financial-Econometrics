%Try one macro varialbe: INFL
INFL_s=(macro(1:(end-1),14)-mean(macro(1:(end-1),14)))/std(macro(1:(end-1),14));
res=ols(er(2:end),[ones(143,1) INFL_s])   %they use standardization of the predictor
res.beta
ans =
    0.0076
    0.0045
res.rsqr
ans =
    0.0127
    
res=ols(er(2:end),[ones(143,1) MS(1:(end-1)) INFL_s])   %they use standardization of the predictor
res.beta
ans =
    0.0076
   -0.0126
    0.0045
res.rsqr
ans =
    0.1100

%A kitchen sink regression does not work (collinearity)
res=ols(er(2:end),[ones(143,1) MS(1:end-1) macro(1:end-1,:)])

%How to get PCs of macro variables (DP, DY, EP,....) in Jiang et al. 2019 JFE data
%Let Y be the columns of macro variables. So Y is 144 by 14

% n=length(Y(:,1));
% Y=macro;
% mu = mean(Y);
% y1 = Y'*Y/n-mu'*mu;  %sample covariance
%[E,v] = eig(y1);    %v is increasing, not what you want
[E,v] = eig(cov(macro));    %v is increasing, not what you want
E'*E  %equals eye(14). I.e. 14 eigen vectors are orthnormal.

[v,ind] = sort(diag(v),'descend');
v  % all eigenvalues

E = E(:,ind(1:5));         %the first five eigenvectors (i.e. PC weights)      %If K=1, I only need to change this line to E = E(:,ind(1:1));  
pcaf = Y*E;    %144 by 5 (5 factors obtained by PCA)


PC1=-(pcaf(:,1)-mean(pcaf(:,1)))/std(pcaf(:,1));
res=ols(er(2:end),[ones(143,1) MS(1:end-1) PC1(1:end-1)])
res.beta
ans =
    0.0076
   -0.0129
    0.0030
res.rsqr
ans =
    0.1030
    
%identify a redundant variable (LTY-TBL=TMS)
[macro(:,9)-macro(:,8),macro(:,11)]

%now only use 12 variables
macro_good=macro(:,[1:3,5:8,10:14]);
res=ols(er(2:end),[ones(143,1) MS(1:end-1) macro_good(1:end-1,:)]);

[E,v] = eig(cov(macro_good));
[v,ind] = sort(diag(v),'descend');
E = E(:,ind(1:12)); 

for j=1:12;
    pcaf = macro_good*E(:,1:j);
    res_pc=ols(er(2:end),[ones(143,1) MS(1:end-1) pcaf(1:end-1,:)]);
    b2(j)=res_pc.beta(2);
end;