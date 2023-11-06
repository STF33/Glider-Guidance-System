%     mdays = ndays_month(dnmb);
%
% Returns # of days in a month/year, diven mtb day or month
% 
function mdays = ndays_month(dnmb);

dv  = datevec(dnmb);
yr  = dv(1);
imo = dv(2);
dn1 = datenum(yr,imo,1);
dn2 = dn1+35;
dv2 = datevec(dn2);
dn2 = datenum(dv2(1),dv2(2),1);

mdays = dn2-dn1;


return
