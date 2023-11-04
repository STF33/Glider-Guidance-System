% Combine RMSE fields for OSE f/casts
% 2 groups: PIES & noPIES
function RMSE = sub_combine_RMSE_OSEfcst(pthout,Z0);
ifc=0;
FCST = {'PIES';'noPIES'};
Nfgr = length(FCST);

for ifc=1:Nfgr
  esimH = FCST{ifc};


  RMSE(ifc).Hnd_name=esimH;
  RMSE(ifc).RMSE_mean=[];
  RMSE(ifc).RMSEprst_mean=[]; % persist

  irun = 0;
  for iyr = 2009:2010
    im1 = 1;
    if iyr == 2009; im1=5; end;
    for imo = im1:12
      frmseout = sprintf('%sRMSE%4.4im_OSEfcst_%s_%i%2.2i.mat',...
                     pthout,abs(Z0),esimH,iyr,imo);
      clear RMSERR
      fprintf('Loading %s\n',frmseout);
      load(frmseout);

      irun = irun+1;
      dmm=RMSERR.ERR_squared;
      dmm=sqrt(nanmean(dmm)); % spatial average
      dmm=dmm(:);
      RMSE(ifc).Fcst_name=esimH;
      ll = size(RMSE(ifc).RMSE_mean,1);
      lld = length(dmm);
      if lld < ll
        dmm(lld+1:ll)=nan;
      elseif ll < lld
        RMSE(ifc).RMSE_mean(ll+1:lld,:) = nan;
      end
      RMSE(ifc).RMSE_mean(:,irun)=dmm;

      dmm=RMSERR.ERRprst_squared;
      dmm=sqrt(nanmean(dmm)); % spatial average
      dmm=dmm(:);
      ll = size(RMSE(ifc).RMSEprst_mean,1);
      lld = length(dmm);
      if lld < ll
        dmm(lld+1:ll)=nan;
      elseif ll < lld
        RMSE(ifc).RMSEprst_mean(ll+1:lld,:) = nan;
      end
      RMSE(ifc).RMSEprst_mean(1:91,irun)=dmm(1:91);
  %   keyboard 
    end
  end
end;


return
