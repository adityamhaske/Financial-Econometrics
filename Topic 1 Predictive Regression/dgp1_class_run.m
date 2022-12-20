%Predictive Regression: Point Estimation

rep=10000;
tic,[b_ols_05,b_bc_05,b_jk2_05,t_ols_05]=dgp1(0,0.5,-0.5,100,rep);toc
tic,[b_ols_09,b_bc_09,b_jk2_09,t_ols_09]=dgp1(0,0.9,-0.5,100,rep);toc
tic,[b_ols_09_m09,b_bc_09_m09,b_jk2_09_m09,t_ols_09_m09]=dgp1(0,0.9,-0.9,100,rep);toc

histogram(b_ols_05,'Normalization','probability')
hold on
histogram(b_ols_09,'Normalization','probability')
histogram(b_ols_09_m09,'Normalization','probability')
