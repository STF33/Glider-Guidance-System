function Tmint = sub_time2mint(TM);
% Convert time mat format
% to:
% minutes since 2009-04-21 12:00:00"
tm0=datenum(2009,4,1,12,0,0);
ddays=TM-tm0;
Tmint = ddays*24*60;


return
