% Combine RMSE fields for forecast group
function RMSE = sub_combine_RMSE_OSSEfcst(IFCST,pthmat,EXPT,Z0);
Nfgr=length(IFCST); % # of forecasts groups 
Ntime = 2;
ifc=0;
for ii=1:Nfgr
  iFcst=IFCST(ii);

  if iFcst >=10; continue; end;

  ifc=ifc+1;
  FCST = sub_fcst_info(iFcst);
  Nruns=1;
  nruns=Nruns;
  Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
  nmexp = FCST.Hindcast_Name;
  pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
  irun1 = FCST.run1;
  irun2 = FCST.run2;
%  ntime = FCST.ntime;
%  ntime = 1;
%  Ntimes = ntime;

  RMSE(ifc).Hnd_name=EXPT(iFcst).Name_short;
  RMSE(ifc).Fcst_nmb=iFcst;
  RMSE(ifc).Time(1).RMSE_mean=[];
  RMSE(ifc).Time(2).RMSE_mean=[];
  RMSE(ifc).Time(1).RMSEprst_mean=[]; % persist
  RMSE(ifc).Time(2).RMSEprst_mean=[]; % persist

  for itime=1:Ntime
    ccn=0;
    RMSE(ifc).Time(itime).RMSE_mean=[];
    for irun=irun1:irun2    % forecast runs
      RUN = FCST.TIME0(itime).RUN(irun);
      pthd1 = RUN.pthbin;
      TM    = RUN.TM;
      YDAY  = RUN.jday;
      nrc   = length(TM);
      DV    = datevec(TM);

      fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
   %    fprintf(' Input data: %s\n',pthd1);

      fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
      frmseout = sprintf('%sRMSE_%s.mat',pthmat,fcstname);
      if (abs(Z0) ~= 10),
        frmseout = sprintf('%sRMSE%4.4im_%s.mat',pthmat,round(abs(Z0)),fcstname);
      end

      clear RMSERR
      fprintf('Loading %s\n',frmseout);
      load(frmseout);

      dmm=RMSERR.ERR_squared;
      dmm=sqrt(nanmean(dmm)); % spatial average
      dmm=dmm(:);
      RMSE(ifc).Time(itime).Fcst_name=fcstname;
      RMSE(ifc).Time(itime).RMSE_mean(:,irun)=dmm;

      dmm=RMSERR.ERRprst_squared;
      dmm=sqrt(nanmean(dmm)); % spatial average
      dmm=dmm(:);
      RMSE(ifc).Time(itime).RMSEprst_mean(:,irun)=dmm;
  %   keyboard 
    end % irun loop - 
  end
end;


return
