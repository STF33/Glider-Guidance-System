% Predictability 1a experiments with NEMO interpolated fields 
% Shifted 1 week each group, each group has +/- 1,2 days wrt to main f/cast
% Grab only main f/cast
function iFcstI = sub_find_predict_fcast(dnmb0,itime);
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
pthmat1 = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';

iFcstI = [];
for ii = 10:16
  FCST = sub_fcstPrdct_info(ii);
  TM = FCST.TIME0(itime).RUN(1).TM;
  if TM(1) == dnmb0;
    iFcstI = ii;
    break;
  end
end

if isempty(iFcstI);
  error('Could not find Predict1a f/cast that starts %s',datestr(dnmb0));
end

return

