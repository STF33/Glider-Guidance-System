function [dS,dE,dC] = get_dtime(dnumb);
%  function [dS,dE,dC] = get_dtime(dnumb);
% Returns day # in NRL hycom convention (NRL):
% First (dS) and last (dE) days of the current month
% and actual day (dC)
% reltaive to day 0 = Dec. 31, 1900
%
% Note that in HYCOM runs, dE = 1st day of the following
% month (to finish the last day), i.e.
% in the limits.* file, dE would be dE+1
% dnumb - matlab date
dH1 = datenum(1900,12,31);
DV  = datevec(dnumb);
dnmS = datenum(DV(1),DV(2),1);
dnmE = dnmS+32;
DV = datevec(dnmE);
dnmE = datenum(DV(1),DV(2),1)-1;

dS = dnmS - dH1;
dE = dnmE - dH1;
dC = dnumb- dH1;

return
