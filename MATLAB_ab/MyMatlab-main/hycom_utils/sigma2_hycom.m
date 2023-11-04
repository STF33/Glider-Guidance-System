       function sgm2 = sigma2_hycom(S,T);
%
% ===============================================
% Caclulate sigma-2 (2000 dB reference pressure)
% as a function of T (deg C) and S (psu)
% based on HYCOM (v. 2.2.27) formula which uses 9-term polynomical fir
% to T [-2:30], S [18:38], modified version of
% Brydon and Sun fit used in HYCOM 2.2.22 version
% see stmt_fns.h
% ===============================================
c1 =  9.903308;     % const. coefficient
c2 = -1.618075e-2;  % T   coefficient
c3 =  7.819166e-1;  % S coef.
c4 = -6.593939e-3;  % T^2 coef.
c5 = -2.896464e-3;  % T S coef.
c6 =  3.038697e-5;  % T^3 coef. 
rc6=  1/c6;
c7 =  3.266933e-5;  % T^2S coef.
c8 =  1.180109e-4;  % S^2 coef.
c9 =  3.399511e-6;  % T S^2 coef. 

pref = 2000.e4;     % Reference pressure;


sgm2 = (c1+S.*(c3+S.*c8)+...
           T.*(c2+S.*(c5+S.*c9)+T.*(c4+S.*c7+T.*c6)));






