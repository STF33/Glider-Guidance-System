% Plot hindcasts (assimilated and freeruns) 
% mean Yucatan vertical section of U field
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

FPLT=[0,0,1,0,0,0,0,0]; % experiments to plot

pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthmatGL = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthtopo_tsis = '/home/ddmitry/codes/HYCOM_TSIS/';
pthtopo = pthtopo_tsis;

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

  if strncmp(esim,esim2,5) | strncmp(esim,esim8,5)
    Tr = sub_avrgUreanls(pthout,pthtopo,esim,ys,ye);
  else
    Tr = sub_avrgUtsis(pthout,pthtopo,esim,ys,ye);
  end

  btx = 'plot_UavYuc_xsct.m';
  fgn=i1;
  kb=3;   % side of the box around Yucatan to plot
  clim1 = -0.3;
  clim2 = 1.2;
  xlim0 = 265;
  ylim0 = -2200;
  dcc = 0.2;
  ctl=sprintf('%s Yuc, Side=%i, %i-%i, dltUcntr=%2.2f',esim,kb,ys,ye,dcc);
  sub_plotUsection_Yucatan(Tr,fgn,kb,ctl,xlim0,ylim0,clim1,clim2,dcc);
  bottom_text(btx,'pwd',1);

end


