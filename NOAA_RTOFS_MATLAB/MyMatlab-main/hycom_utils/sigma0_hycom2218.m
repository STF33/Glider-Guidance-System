       function sgm0 = sigma0_hycom2218(S,T);
%
% ===============================================
% Caclulate sigma-0 (0 dB reference pressure)
% as a function of T (deg C) and S (psu)
% based on HYCOM (v. 2.2.18) formula which uses  polynomical fit
% to T [-2:30], S [18:38]
% see stmt_fns.h
% ===============================================
c1 = -1.36471e-01;  % const. coefficient
c2 =  4.68181e-2;   % T   coefficient
c3 =  8.07004e-1;  % S coef.
c4 = -7.45353e-3;  % T^2 coef.
c5 = -2.94418e-3;  % T S coef.
c6 =  3.43570e-5;  % T^3 coef. 
rc6=  1/c6;
c7 =  3.48658e-5;  % T^2S coef.


pref = 0.0;     % Reference pressure;
sgm0 = c1+c3*S+T.*(c2+c5*S+T.*(c4+c7*S+c6*T));
