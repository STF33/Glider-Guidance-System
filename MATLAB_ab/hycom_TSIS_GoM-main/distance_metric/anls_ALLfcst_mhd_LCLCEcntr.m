% Analyze MHD for all forecasts
% 
% Before analysis:
% extract LC/LCE contour
% compute MHD
% 
% Calculated MHD for the LC contours in 
% mhd_LCLCEcntr_nemo_fcsthycom.m
% LC is extracted in hycom_TSIS/extr_lc_hycom_nemoV1.m
% NEMO LC is extracted in hycom_TSIS/extr_lc_ssh_nemo.m
%
% Hindcasts used for the forecast runs:
% H/cast #2 - Full 2D SSH                         not finished
%        #3 - AVISO SSH tracks only               Fcst 
%        #6 - 30th pnt T/S profiles GoM           not finished
%        #7 - AVISO + UGOS PIES (T/S profiles)    Fcst 
%        #8 - AVISO + extended PIES               Fcst 
% 
% Plot time series:
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;

startup;

close all
clear

% Figure to plot:
f_mhd_tser = 0; % plot MHD time series 


pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst/';
btx = 'anls_ALLfcst_mhd_LCLCEcntr.m';

IFCST=[2,3,6,7,8]; % Frecast groups to analyze
Nfgr=length(IFCST); % # of forecasts groups 

% Hindcast info
load('hycom_tsis_expts.mat');  % EXPT array

ixx=0;
ymx = 0;
for ii=1:Nfgr
  iFcst=IFCST(ii);

  FCST = sub_fcst_info(iFcst);
  Nhind = FCST.Nhind;  % initial cond from the  hindcast #
  hnd_name = FCST.Hindcast_Name;
  pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields
  irun1 = FCST.run1;
  irun2 = FCST.run2;
  ntime=FCST.ntime; % 2 time windows for forecasts
  Ntimes=ntime;

  for itime=1:ntime  % 2011 and 2012 time windows
    for irun=irun1:irun2
      nmexp   = sprintf('fcst%2.2i-%2.2i%2.2i',Nhind,itime,irun);
      fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycom%s.mat',...
                  pthmat,nmexp);
      fprintf('Loading %s\n',fmat1);
      A=load(fmat1);  % MHD for f/cast and persistence
      DV=datevec(A.TM);
      ixx=ixx+1;
      MHD(ixx).mhd = A.MHD(:,1);
      MHD(ixx).Name = nmexp;
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).nemo0 = A.MHD(:,2); % NEMO true contour on day 1
      MHD(ixx).hycom0= A.MHD(:,3); % HYCOM contour on day 1
      MHD(ixx).iFcst=iFcst;    % Forecast group # (groupped by hindcast initial fields)
      MHD(ixx).Hindcast_Number=Nhind;
      MHD(ixx).Time_period=itime;
      MHD(ixx).Run_number=irun;
    end
  end
end
Ntot = ixx;

% Find corresponding hindcasts:
for ifc=1:Nfgr
  for ii=1:length(MHD);
      if MHD(ii).iFcst==IFCST(ifc); break; end;
  end
  HDCST(ifc)=MHD(ii).Hindcast_Number;
end

% Hindcast colors 
% match with hindcast for
% comparison
CLR = [0.4 0.8 0.2;...
       0   0.6 0.9; ...
       0.5   1   0.7; ...
       1   0.4 0.5; ... 
       0.  1  0.3; ...
       1   0   0.9; ...
       0.8 0   0.4; ...
       1   0.8 0; ...
       0.8 0.5 0; ...
       0.7 0.6 0.4; ...
       0.5 0.2 1];

% Persistence:
%ipst = Ntot+1;
%CLR(ipst,:) = [0.4 0.4 0.4];


%MXD = A.MAXD;
TM=A.TM;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);


%keyboard
% ----------------------
%
%  Plot time series of MHD
%
% -----------------------
ifn = 0;
clear ixx

if f_mhd_tser==1
  for nhnd=7:7
    for itime=1:Ntimes  % 2011 and 2012 time windows
      for irun=1:7
        ifn=ifn+1;
        itot = (nhnd-7+itime-1)*Nruns+irun;

        sub_plot_mhd_tser(MHD,CLR,irun,ifn,ipst,itot);
        bottom_text(btx,'pwd',1,'Position',[0.09 0.02 0.4 0.04]);

  %
  % Calculate MHD cum score
  %      mhd=MHD(itot).mhd;
  %      mhdP=MHD(itot).hycom0;
  %      MHD(itot).cumMHD=cumsum(mhd);
  %      MHD(itot).cumPST=cumsum(mhdP);
  %     keyboard
      end  % runs
    end
  end
end

% Weekly pools
WK=[1:7:92]; % weekly pools
nwk=length(WK)-1;

% Mean MHD by forecast-hindcast/time groups
nll = length(MHD);
cMHD=[];
cMHDp=[];
irun=0;
iFcst=[];
for ifc=1:Nfgr
  for itime=1:2
    POOL(ifc).Time(itime).pm1=[];  % month 1 pool
    POOL(ifc).Time(itime).pm2=[];
    POOL(ifc).Time(itime).pm3=[];

    for iwk=1:nwk
      POOL(ifc).Time(itime).WK(iwk).pw=[];
    end
  end
end


for ill=1:nll
  if isempty(iFcst); iFcst=MHD(ill).iFcst; end
  if MHD(ill).iFcst ~= iFcst;
    irun=0;
    iFcst=MHD(ill).iFcst;
  end
  irun=irun+1;
  ifc=find(IFCST==iFcst);
  itime=MHD(ill).Time_period;
% HYCOM mean MHD by 30 days
  dmm=MHD(ill).mhd;
%  cMHD(ifc,irun,1)=nanmean(dmm(1:30));
%  cMHD(ifc,irun,2)=nanmean(dmm(31:60));
%  cMHD(ifc,irun,3)=nanmean(dmm(61:end));
%
% Pool all runs for the same forecast group and time period
  pm1=POOL(ifc).Time(itime).pm1;
  pm2=POOL(ifc).Time(itime).pm2;
  pm3=POOL(ifc).Time(itime).pm3;

% Monthly pools
  pm1=[pm1;dmm(1:30)];
  pm2=[pm2;dmm(31:60)];
  pm3=[pm3;dmm(61:end)];

  POOL(ifc).Time(itime).pm1=pm1;
  POOL(ifc).Time(itime).pm2=pm2;
  POOL(ifc).Time(itime).pm3=pm3;

% Weekly pools
  for iwk=1:nwk
    id1=WK(iwk);
    id2=WK(iwk+1)-1;

    pw=POOL(ifc).Time(itime).WK(iwk).pw;
    pw=[pw;dmm(id1:id2)];
    POOL(ifc).Time(itime).WK(iwk).pw=pw;
  end  

end


%keyboard

POOL(1).ylim=[0,160];
% ----------------------------
% Plot: all F/csts in 1 group for 3 30-day time periods
% coordinates of the Bars relative to X=# of the 30-day period, 1,2,3
nfg=40;
anls_nm='MHD(km) HYCOM LCLCE months';
sub_plotMHD_bars_month(nfg,Nfgr,IFCST,CLR,POOL,EXPT,anls_nm);
bottom_text(btx,'pwd',1);  

% Plot weekly bars
nfg=41;
anls_nm='MHD(km) HYCOM LCLCE weeks';
sub_plotMHD_bars_week(nfg,Nfgr,IFCST,CLR,POOL,EXPT,WK,anls_nm);
bottom_text(btx,'pwd',1);





keyboard

%
% Plot skill vs initial error  
% Group by forecast groups
SKL=[];
ntot=length(MHD);
for jj=1:Nfgr
  for itime=1:2
    SKL(jj).Time(itime).MHD=[];
  end
end

for itot=1:ntot
  iFcst=MHD(itot).iFcst;
  jj=find(IFCST==iFcst);
  itime=MHD(itot).Time_period;
  irun=MHD(itot).Run_number;
  mhd=MHD(itot).mhd;
  SKL(jj).Time(itime).MHD(irun,1)=mhd(1);
  SKL(jj).Time(itime).MHD(irun,2)=median(mhd(2:31));
  SKL(jj).Time(itime).MHD(irun,3)=median(mhd(32:61));
  SKL(jj).Time(itime).MHD(irun,4)=median(mhd(62:91));
  amm(itot,1)=mhd(1);
%  SKL(jj).Time(itime).MHD(irun,2)=mhd(2:31);
end


dxp=0.25;
POST=[0.3 0.67 dxp dxp; ...
      0.3 0.37 dxp dxp; ...
      0.3 0.09 dxp dxp];
 
for itime=1:2
  figure(40+itime); clf;

  for imean=1:3 ;  % 1-30, 31-60, 61-90-day mean skills
    pss=POST(imean,:);
    axes('Position',pss);
    hold on;

    mxy=0;
    mxx=0;
    mnx=1e6;
    for jj=1:Nfgr  % forecasts
      jjr=imean+1; 
      iFcst=IFCST(jj);
      clr=CLR(iFcst,:);

      skl=SKL(jj).Time(itime).MHD;
      if ~isempty(skl)
      plot(skl(:,1),skl(:,jjr),'.','Markersize',18,'Color',clr);
      mxy=max([mxy,max(skl(:,jjr))]);
      mnx=min([mnx,min(skl(:,jjr))]);
      end
    end
    mxx=max(skl(:,1));
    xl2=1.1*max([mxx,mxy]);
    xl1=0.8*min([mnx,min(skl(:,1))]);
   
    axis('equal'); 
    set(gca,'tickdir','out',...
            'xlim',[0 xl2],...
            'ylim',[0 xl2],...
            'Fontsize',16,...
            'xgrid','on',...
            'ygrid','on');

    if itime==1,
      str='01/05/2011-12/06/2011';
    else
      str='01/01/2012-12/02/2012';
    end
    if imean==1
      dtr='1-31';
    elseif imean==2
      dtr='32-61';
    else
      dtr='62-91';
    end
    spp=sprintf('Median(MHD) vs MHD(t=0), %s, %s',str,dtr);
    title(spp);
  end
end







