function dnmb = inverse_hycom_dtime(hnmb);
% function dnmb = inverse_hycom_dtime(hnmb);
% inverse hycom date to normal date in maltab convention
%
% dnumb - matlab date
dH1 = datenum(1900,12,31);
dnmb = dH1+hnmb;


return
