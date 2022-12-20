function [b_hjal,white_t,hodrick_FR,hodrick_mod]=hodrick_tstat_10_M524(r,x,h,R_null,b_null);

%from hodrick_tstat_10.m

%x is n by K
%R_null: q by K
%b_null: q by 1
%the null is R_null*b=b_null
%marginal predictabitliy: R_null=[0 1];b_null=0;

%e.g. [b_hjal,white_t,hodrick_FR,hodrick_mod]=hodrick_tstat_10_M524([mean(ret_4);ret_4],dp_4,3,1,0), sqrt(hodrick_FR)

%Output:

%b_hjal: point estimate (horizon-h slope in the forward regression)
%hodrick_FR: t-stat of horizon-h forward regression, using Hodrick SE (known as Hodrick 1B standard error)

n=length(r);
K=length(x(1,:));

iota=ones(n-1,1);  %regression of r
res1=ols(r(2:n),[iota x(1:(n-1),:)]);
b0=res1.beta(1);
b=res1.beta(2:end);
u0=res1.resid;

%Part 1: t-test with Hodrick 1B standard error
y=zeros(n-h+1,1);
for i=1:(n-h+1);
    s=0;
    for j=0:(h-1);s=s+r(i+j);end;
    y(i)=s;
end;

iota=ones(n-h,1);
res=ols(y(2:(n-h+1)),[iota,x(1:(n-h),:)]);
b_hjal=res.beta(2:end);

%Hodrick stat
wk=zeros(n-h,K+1);
XX=[ones(n,1) x];
for i=1:(n-h);
    gj=0;
    for j=1:h;gj=gj+XX(i+h-j,:);end;
    wk(i,:)=u0(i+h-1)*gj;
end;

meat=(wk'*wk);
bread=(XX(1:(n-h),:)'*XX(1:(n-h),:))^(-1);

se_hod=bread*meat*bread;
se_hod=se_hod(2:end,2:end);

hodrick_FR=(R_null*b_hjal-b_null)'*(R_null*se_hod*R_null')^(-1)*(R_null*b_hjal-b_null);

%Modified Hodrick stat
hodrick_mod=0;

%White's t-stat, ignoring serial correlation
iota=ones(n-h,1);
res=hwhite(y(2:(n-h+1)),[iota,x(1:(n-h),:)]);
white_t=res.tstat(2:end);