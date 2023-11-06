       function sgm2 = sigma2_hycom(S,T);
%
% ===============================================
% Caclulate sigma-2 (2000 dB reference pressure)
% as a function of T (deg C) and S (psu)
% based on HYCOM formula which uses 7-term polynomical fit
% to T [-2:30], S [18:38], 
% Brydon and Sun fit (HYCOM 2.2.22 version)
% see stmt_fns.h
% ===============================================
c1 =  9.77093;     % const. coefficient
c2 = -2.26493e-2;  % T   coefficient
c3 =  7.89879e-1;  % S coef.
c4 = -6.43205e-3;  % T^2 coef.
c5 = -2.62983e-3;  % T S coef.
c6 =  2.75835e-5;  % T^3 coef. 
c7 =  3.15235e-5;  % T^2S coef.

pref = 2000.e4;     % Reference pressure;


sgm2 = c1+S.*c3+...
           T.*(c2+c5*S+T.*(c4+S.*c7+T.*c6));






