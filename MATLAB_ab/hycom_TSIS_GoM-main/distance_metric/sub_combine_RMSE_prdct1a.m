% Combine RMSE from Predictability 1a run 1 (not shift of IC +/- 2,1 days)
% Update RMSE for OSSE f/casts 2,3,6,7,8 by adding f/cast 10,11,..., 16 - interpolated NEMO
% Predictability 1a forecasts - initialized from NEMO+GLORYS 
% interpolated onto HYCOM horiz/vert grid 
% 7 forecasts ran with 1 week shift starting May 1
% combine all 7 runs as 1 forecast group to compare 
% with OSSE forecast groups (AVISO, PIES, etc.)
% mhd - see mhd_LCLCEcntrPrdct_nemo_fcsthycom.m
% Grab only runs 1 e.g. MHD_LCLCE_nemo_persist_OSSEfcst16-0201.mat
% 0202 is IC with 2day shit, 0203 -1 day, 0204 +1, 0205 +2 - not needed

function RMSE = sub_combine_RMSE_prdct1a(RMSE,pthout,Z0); 
IFCST = [10:16];
ifc   = length(RMSE);
ifc   = ifc+1;
Ntime = 2;
Nfgr  = length(IFCST); % # of forecasts groups 
lng   = 91; 

for itime = 1:2
  irr = 0;
  for iFcst = 10:16
    FCST  = sub_fcstPrdct_info(iFcst);
    Nhnd  = FCST.Nhind;  % initial cond from the  hindcast #
    nmexp = FCST.Hindcast_Name;
    pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields

    irun = 1;  % pick only prdct 1a runs
 
    fcstname = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
    frmseout = sprintf('%sRMSE_%s.mat',pthout,fcstname);
    if (abs(Z0) ~= 10),
      frmseout = sprintf('%sRMSE%4.4im_%s.mat',pthout,round(abs(Z0)),fcstname);
    end

    clear RMSERR
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    irr = irr+1;
    dmm=RMSERR.ERR_squared;
    dmm=sqrt(nanmean(dmm)); % spatial average
    dmm=dmm(:);
    RMSE(ifc).Time(itime).Fcst_name=fcstname;
    RMSE(ifc).Time(itime).RMSE_mean(1:lng,irr)=dmm(1:lng);

    dmm=RMSERR.ERRprst_squared;
    dmm=sqrt(nanmean(dmm)); % spatial average
    dmm=dmm(:);
    RMSE(ifc).Time(itime).RMSEprst_mean(1:lng,irr)=dmm(1:lng);
  end
end;


return

