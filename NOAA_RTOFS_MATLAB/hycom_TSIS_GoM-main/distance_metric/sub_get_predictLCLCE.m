% Load mat files with LC/LCEs contours 
% for predict 1a
function [LCE,LCXY] = sub_get_predictLCLCE(iFcst,itime,irun);

pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
FCST = sub_fcstPrdct_info(iFcst);

flnm = sprintf('hycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',iFcst,itime,irun);
fmatout = sprintf('%s%s',pthmat2,flnm);

fprintf('Loading %s\n',fmatout);
load(fmatout);

TM  = LCXY.TM;
fprintf('Predict1a day1=%s\n',datestr(TM(1)));


return
