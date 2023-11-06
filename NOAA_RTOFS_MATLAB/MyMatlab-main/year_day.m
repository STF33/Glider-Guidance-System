     function yrday = year_day(dnum);
%     function yrday = year_day(dnum);
%
% Calculates year day of the day given in matlab day number dnum
%

M     = datevec(dnum);
dJ1   = datenum(M(1),1,1); 
yrday = dnum-dJ1+1;


