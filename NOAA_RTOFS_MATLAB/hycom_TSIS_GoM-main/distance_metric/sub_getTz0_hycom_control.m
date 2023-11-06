% Get interpolated onto Z0 T field
% from the control run
function FF = sub_getTz0_hycom_control(dnmb,iFcst,itime,irun,Z0);

pth1 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
flnm = sprintf('%shycom_t2Z%4.4i_fcstPrdct%2.2i-%2.2i%2.2i.mat',...
                pth1,abs(Z0),iFcst,itime,irun);

fprintf('Loading control run: %s\n',flnm);
load(flnm);

%keyboard

TM = TZH.TM;
iTm = find(TM==dnmb);
if isempty(iTm),
  error('Control run no date %s',datestr(dnmb));
end

% HYCOM field:
FF = squeeze(TZH.Tz(iTm,:,:));

return
