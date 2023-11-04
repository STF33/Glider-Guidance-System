% Calculate nest date for given year, month, skipping # of days
function [dS, dE] = sub_nest_dates(yr,mo,mday,dday);
dH1 = datenum(1901,1,1);
dm1 = datenum(yr,mo,mday)+dday-1;
ndays = [dH1:dday:dm1];
dS = ndays(end-1);
dE = ndays(end);
fprintf('Nest days for requested %s are %s & %s\n',datestr(datenum(yr,mo,mday)),...
				datestr(dS), datestr(dE));

return



