% Compute RMSE for Yuc transport wrt to GLBu run
% that is used as nest fields for HYCOM-TSIS 
%
% calculated in
% calc_volFlux_GLBb.m - Global reanalysis
% and calc_volFlux_v3.m for HYCOM-TSIS
% and calc_volFlux_GoMreanalysis.m for 3hr output 0.04 GoM reanls with tides
%
% Transport saved by years
%
% In order to match GLBu reanalysis Yucatan Transport
% the following adjustment have been applied:
% Barotropic transport increased by 2.5 times
% Montgomery potential is regressed using GLB reanalysis fields
% mean dynamic topography in the data/observation files (used for assimilation)
% is replaced with the mean SSH map calculated from the
% GLB reanalysis
% 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all

esim1='PIES'; % adjusted nest&Montg hindcast data assimilative runs with deep PIES obs
esim2='GLBu191'; % 
esim3='noPIES'; % adjusted nest hindcast data assimilative runs without deep PIES obs
esim4='freerun_noadj'; % freerun with original nest fields
esim5='PIESnoadj'; % original hindcast data assimilative runs with deep PIES obs
esim6='noPIESnoadj'; % original hindcast data assimilative runs without deep PIES obs
esim7='GLfreerun';  % GLORYS freerun
esim8='GoMu501'; % 0.04 NCODA + GoMu reanalysis

FPLT=[1,1,1,0,0,0,0,1]; % experiments to plot
iG = 2; % GLB run

pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmatGL = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

xsct_name = 'Yucatan';

ys=2009;
ye=2011;

cc=0;

II = find(FPLT==1);
for i1=1:length(II)
  jj=II(i1);
  switch(jj)
   case(1)
    esim=esim1;
    pthout = pthmat;
   case(2)
    esim=esim2;
    pthout = pthmat;
   case(3)
    esim=esim3;
    pthout = pthmat;
   case(4)
    esim=esim4;
    pthout = pthmat;
   case(5)
    esim=esim5;
    pthout = pthmat;
   case(6)
    esim=esim6;
    pthout = pthmat;
   case(7)
    esim=esim7;
    pthout = pthmatGL;
   case(8)
    esim=esim8;
    pthout = pthmat;
  end

  cc=cc+1;

  Tr=sub_combine_transp(pthout,esim,ys,ye);
  TR(cc).Name=esim8;
  TR(cc).YR=Tr.YR;
  TR(cc).Tb=Tr.Barotrop;
  TR(cc).Ttot=Tr.Total;
  TR(cc).Time=Tr.Time;

end

nTR=cc;

TM0=TR(iG).Time;
TV0=TR(iG).Ttot;

for i1=1:length(II)
  jj=II(i1);
  if jj==iG; continue; end;

  TM=TR(i1).Time;
  TV=TR(i1).Ttot;

  nTM = length(TM);
  cc = 0;
  rsum = 0;
  for it=1:nTM
    ix = find(round(TM0)==round(TM(it)));
    if isempty(ix), continue; end;
    if isnan(TV(it)), continue; end
    cc=cc+1;
    rsum = rsum+(TV0(ix)-TV(it))^2;
  end

  RMSE(i1)=sqrt(rsum/cc);

end









