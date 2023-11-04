% Combine RMSE fields for OSE hindcasts
% 2 groups: PIES & noPIES vs GOMu
% for 2009 and 2011 - not full year - only covered by PIES obs

function RMSE = sub_combine_RMSE_OSEhndcst(pthout,Z0);
ifc=0;
FCST = {'PIES';'noPIES'};
Nfgr = length(FCST);

for ifc=1:Nfgr
  esimH = FCST{ifc};


  RMSE(ifc).Hnd_name=esimH;
  RMSE(ifc).RMSE_mean=[];
  RMSE(ifc).RMSEprst_mean=[]; % persist

  icc = 0;
  for iyr = 2009:2011
    frmseout = sprintf('%sRMSE%4.4im_OSEhindcast%s_GOMu_%i.mat',...
                   pthout,abs(Z0),esimH,iyr);
    clear RMSERR
    fprintf('Loading %s\n',frmseout);
    load(frmseout);

    icc = icc+1;

    dmm=RMSERR.ERR_squared;
    dmm=sqrt(nanmean(dmm)); % spatial average
    dmm=dmm(:);
    if iyr == 2009
      dmm = dmm(91:end);  % Apr - Dec
    elseif iyr == 2011
      dmm = dmm(1:314);  % Jan - Nov.
    end
    RMSE(ifc).Fcst_name=esimH;
    RMSE(ifc).Year(icc) = iyr;
%
% Annual median, IQR:
    mdn = median(dmm);
    lup = prctile(dmm,75);
    lbt = prctile(dmm,25);   
    RMSE(ifc).median(icc) = mdn;
    RMSE(ifc).p25(icc)    = lup;
    RMSE(ifc).p75(icc)    = lbt;
  end
end;


return
