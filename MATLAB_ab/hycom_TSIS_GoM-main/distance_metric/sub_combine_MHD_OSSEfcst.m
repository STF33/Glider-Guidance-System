% Combine MHD scores for OSSE forecasts
% These are experiments # < 10
function MHD = sub_combine_MHD_OSSEfcst(IFCST,pthmat);

Nfgr=length(IFCST); % # of forecasts groups
ifc = 0;
Ntime = 2;
for ii=1:Nfgr
  iFcst=IFCST(ii);

  if iFcst >=10; continue; end;

  ifc=ifc+1;
  FCST = sub_fcst_info(iFcst);
  Nhnd = FCST.Nhind;  % initial cond from the  hindcast #
  nmexp = FCST.Hindcast_Name;
  pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
  irun1 = FCST.run1;
  irun2 = FCST.run2;
  ntime = FCST.ntime;


  MHD(ifc).Fcst_OSSE = nmexp;
  MHD(ifc).Fcst_nmb  = iFcst;
  MHD(ifc).Time(1).MHD_mean=[];
  MHD(ifc).Time(2).MHD_mean=[];

  for itime=1:Ntime
    MHD(ifc).Time(itime).MHD_mean=[];
    for irun=irun1:irun2    % forecast runs
      RUN = FCST.TIME0(itime).RUN(irun);
      pthd1 = RUN.pthbin;
      TM    = RUN.TM;
      YDAY  = RUN.jday;
      nrc   = length(TM);
      DV    = datevec(TM);

      fcst_name = sprintf('%s%2.2i-%2.2i%2.2i',nmexp,iFcst,itime,irun);
      fprintf(' Forecast: %i, TimePeriod %i, FcstRun %i\n',Nhnd,itime,irun);
      fmhdout = sprintf('%sMHD_LCLCE_nemo_persist_OSSEfcst%2.2i-%2.2i%2.2i.mat',...
                  pthmat,Nhnd,itime,irun);

      fprintf('Loading %s\n',fmhdout);
      A=load(fmhdout);

      dmm = A.MHD;
      MHD(ifc).Time(itime).Fcst_name   = fcst_name;
      MHD(ifc).Time(itime).MHD(:,irun) = dmm(:,1);
      MHD(ifc).Time(itime).MHD_prst(:,irun) = dmm(:,2);

    end % irun loop
  end
end


return
