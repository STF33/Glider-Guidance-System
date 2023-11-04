% Read pies T/S profiles
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

YR=2011;

% Synthetis OSSE observations
% frmo NEMO 1km free run 2011
pthosse = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/OSSE/';
pthout  = sprintf('%spies/',pthosse);
pthpies = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/pies/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

% Read the topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


%fnm = sprintf('%sia511.nc',pthpies);

%
% Create OSSE T/S netCDF files
finp=sprintf('%sugos_mooring_%i.mat',pthosse,YR);
fprintf('Loading %s\n',finp);
A=load(finp);
plat=A.ugos_lat;
plon=A.ugos_lon;
pZ=double(A.ugos_depth);
np=length(plat);

IJ = sub_XY2indx([plon,plat],LON,LAT);


% Extended moorings:
finpE=sprintf('%sextd_mooring_2011.mat',pthosse);
A2=load(finpE);
platE=A2.extd_lat;
plonE=A2.extd_lon;
np2=length(platE);

f_plt=0;
if f_plt==1
		figure(1); clf;
		contour(LON,LAT,HH,[0 0],'k');
		hold on;
		contour(LON,LAT,HH,[-5000:1000:-10],'Color',[0.6 0.6 0.6]);
		axis('equal');
  plot(plonE,platE,'+','Color',[0 0.5 0.8]);
		plot(plon,plat,'r.','Markersize',16);
		set(gca,'tickdir','out',...
										'xlim',[-96 -80],...
										'ylim',[19 30.5]);
		title('UGOS & Extended, OSSE PIES locations');
		btx='create_pies_TS.m';
		bottom_text(btx,'pwd',1);
end

npies=length(plat);

% Extend bottom
for iip=1:npies
  flnc=sprintf('%sugos%3.3i.nc',pthout,iip);
  TT=squeeze(A.ugos_temp(iip,:,:));
  SS=squeeze(A.ugos_salt(iip,:,:));
  
  ib=min(find(TT(:,1)==0));
  TT(ib,:)=TT(ib-1,:);
  SS(ib,:)=SS(ib-1,:);

  A.ugos_temp(iip,:,:)=TT;
  A.ugos_salt(iip,:,:)=SS; 
 
end




TM=A.ugos_time;
TM=TM(:);
Tmint = sub_time2mint(TM);
ntime=length(Tmint);

for iip=1:npies
  fprintf('Location %i\n',iip);
  fnmb=sprintf('ugos%3.3i',iip);
%  fcdl=fprintf('%sugos%3.3i.cdl',pthout,fnmb);
  sub_cdl(fnmb,pthout,ntime);

  flnc=sprintf('%s%s.nc',pthout,fnmb);

  T=squeeze(A.ugos_temp(iip,:,:));
  S=squeeze(A.ugos_salt(iip,:,:));

  nc_varput(flnc,'depth',pZ);
  nc_varput(flnc,'time',Tmint);
  nc_varput(flnc,'T_var',T');
  nc_varput(flnc,'S_var',S');
end


% Check the bottom:
for iip=1:npies
  flnc=sprintf('%sugos%3.3i.nc',pthout,iip);
  TT=nc_varget(flnc,'T_var');
  TT=TT';

  ib=min(find(TT(:,1)==0));
  hb=-pZ(ib-1);

  hbH=HH(IJ(iip,2),IJ(iip,1));
  if abs(hb)<abs(hbH)
    fprintf('Above HYCOM bottom: %i hb=%6.1f HYCOM=%6.1f\n',iip,hb,hbH);
  end

end



f_chck=0;
if f_chck==1
  iip=10;
  flnc=sprintf('%sugos%3.3i.nc',pthout,iip);
  TT=nc_varget(flnc,'T_var');
  TT=TT';

  T=squeeze(A.ugos_temp(iip,:,:));
  Z=-pZ;   

  it=20;
  figure(1); clf;
  plot(TT(:,it),Z,'.-');
  hold on;

  hb=HH(IJ(iip,2),IJ(iip,1));
  plot([3 6],[hb hb],'b-');
  

end





