function mntg1 = sub_montg1_regr(ssh);
% Use regression to estimate montg1
% from the SSH field, m
% see: fit_regr_montgomery.m

gg     = 9.806;
pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/080/data_anls/');
fmat=sprintf('%sregression_montg_ssh.mat',pthmat);

load(fmat);

B0=RGR.intercept;
B1=RGR.slope;

mntg1 = B0 + B1.*ssh;

return