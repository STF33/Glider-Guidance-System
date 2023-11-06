% GLORYS analysis 
% 6-day forecasts
function [flnm,flgr] = sub_find_GLORYS_file(pthglorys,dnmb);
dE1 = datenum(2011,1,5);

DE = [dE1:7:datenum(2012,12,31)]; % end of forecast windows
iM = min(find(DE>dnmb));

dv1 = datevec(dnmb);
dv2 = datevec(DE(iM));
flnm = sprintf('mercatorglorys12v1_gl12_mean_%4.4i%2.2i%2.2i_R%4.4i%2.2i%2.2i.nc',...
                dv1(1:3),dv2(1:3));

flgr = sprintf('%s%s',pthglorys,flnm);

ie = exist(flgr,'file');
if ie<=0
  fprintf('Does not exist %s\n',ie);
end


return
