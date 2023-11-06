% Plot SSH from OSSE hindcasts 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hausdorff_dist;
startup;

close all
clear

% Set flags for extracting experiments:
EXON = zeros(9,1);
EXON(2) = 1; % select expt to plot
dnmb = datenum(2011,1,5);  % date to plot

huge = 1.e20;
rg=9806;  % convert pressure to depth, m
% Hindcasts:
pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

%fmatout = sprintf('%sLC_coord_osse_hycomV1.mat',pthmat);

btx = 'plot_ssh_OSSEhndcst.m';

fhnd = 'hycom_tsis_expts.mat';
fprintf('Loading %s\n',fhnd);
load(fhnd);

Nruns = length(EXON);

for ii=1:Nruns
  if EXON(ii)==0;
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
end


for ii=1:Nruns
 fprintf('%i: %s \n',ii,EXPT(ii).path);
end


% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


% GoM region HYCOM:
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM


% section to plot:
SCT = sub_GoM_xsct(LON,LAT,HH);



%cmp=flipud(colormap_cold(360));
% Colormap
ncc=100;
cl1=[0,0.3,0.6];
cl2=[1,1,1];
clr1=mix_2colors(cl1,cl2,ncc);

%cl1=cl2;
%cl2=[0.6,1,0.8];
cl1=cl2;
cl2=[0.6,0.9,0.7];
clr2=mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,0.8,0.2];
clr3=mix_2colors(cl1,cl2,ncc);
cmp = [clr1;clr2;clr3];
cmp = smooth_colormap(cmp,5);



dv = datevec(dnmb);
yr = dv(1);
mo = dv(2); 
dm = dv(3);
iday = dnmb-datenum(yr,1,1)+1;
DV = dv;


% PLOTTING
Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;
  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fprintf('%s %s\n',nmexp,fmatout);

	fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
	finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
	fin=fina;

	ie = exist(fin,'file');

	if ~ie 
		fprintf('Missing: %s\n',fin);
		continue;
	end 

	fprintf('Reading %s\n',fina);
	fld = 'srfhgt';
	[F,nn,mm,ll] = read_hycom(fina,finb,fld);
	F(F>huge)=nan;
	ssh=squeeze(F)./(1e-3*rg);  % ssh m
% 
% Subtract anomaly:
	dmm=ssh;
	dmm(INH==0)=nan;
%  dmm(HH>-200)=nan;
	sshM=nanmean(nanmean(dmm));
	ssh=ssh-sshM;


  figure(jj); clf;
  set(gcf,'Position',[1573         560         949         751]);
  axes('Position',[0.08 0.4 0.8 0.5]);
   % hindcast
  pcolor(LON,LAT,ssh); shading flat;
  colormap(cmp);
  clm1=-0.3;
  clm2=0.6;
%  caxis([-0.5 0.5]);
  caxis([clm1, clm2]);
  hold;

  contour(LON,LAT,ssh,[0.17 0.17],'k-','linewidth',1.6);

% Plot xsection line:
  XX=SCT.long;
  YY=SCT.latd;
  plot(XX,YY,'r.');

% Plot Distance hash-lines
	DST  = cumsum(SCT.dX);
	xplt = [0:200:max(DST)+10];
	for idx=1:length(xplt)
		L0=xplt(idx);
		dlt= abs(DST-L0);
		jd = find(dlt==min(dlt));
		xl0= XX(jd);
		yl0= YY(jd);
		if idx==1
			plot(xl0,yl0,'k.','Markersize',20);
		else
			plot(xl0,yl0,'r.','Markersize',20);
		end
		
	end



  axis('equal');
  set(gca,'xlim',[-97.8 -80.7],...
   'ylim',[18. 31],...
   'xtick',[-98:2:-82],...
   'ytick',[18:2:32],...
   'tickdir','out',...
   'color',[0 0 0],...
   'Fontsize',12);
  hb=colorbar;
  set(hb,'Position',[0.81 0.35 0.015 0.5],...
   'Ticks',[-1.0:0.2:1.0],...
   'Fontsize',12);
  
 % stl=sprintf('ssh, NEMO, LC/LCE contours (17cm) HYCOM hindcast %2.2i, %s',ixx-1,datestr(dnmb));
  stl=sprintf('ssh, HYCOM OSSEhindcast %s %s',nmexp,datestr(dnmb));
  title(stl);

  bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.4 0.04]);
%keyboard
end



  
  
  
  
  



