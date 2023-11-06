% Check created nest archive file
% Plot SSH NEMO vs HYCOM expt
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear
close all

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

nlrs = 30;     % depth layers in HYCOM TSIS

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthintrp  = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthmat    = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/datamat/';


iyr = DV(1);
iday = dnmb-datenum(iyr,1,1)+1;

flnm = sprintf('archv.%4.4i_%3.3i_00.a',iyr,iday);
fina = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthoutp,iyr,iday);
finb = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthoutp,iyr,iday);

%
% NEMO:
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);
[mm,nn] = size(LONN);
%[DXN,DYN]=sub_dx_dy(LONN,LATN);
[XM,YM]=meshgrid([1:nn],[1:mm]);
% GoM region, NEMO:
GOMN=[         100         365
         651         337
        1091         687
        1246         798
        1512         881
        1787         998
        1904        1292
        1710        1914
          23        1920
           8         748];
INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
clear XM YM



% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;

%Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mh,nh]=size(HH);
m=mh;
n=nh;
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

dmean = 1;  % = 1 - subtract spatial mean within GoM
sshn = sub_getSSH_nemo(DV,LONN,LATN,INN,dmean);


% HYCOM SSH:
[F,n,m,l] = read_hycom(fina,finb,'srfhgt');
F = squeeze(F);
F(F>huge)=nan;
sshh=squeeze(F)./(1e-3*rg);  % ssh m

% Subtract anomaly:
dmm=sshh;
dmm(INH==0)=nan;
%  dmm(HH>-200)=nan;
sshM=nanmean(nanmean(dmm));
sshh=sshh-sshM;

% Derive LC contour:
% 
dmm=sshh;
dmm(INH==0)=nan;
%  LCH = identify_LC(LON,LAT,dmm,Bisol);


fprintf('Plotting ...\n');
% =====================
% Plotting
% =====================
%  cntr=cntr+1;
% NEMO:
yr = DV(1);
mo = DV(2);
dm = DV(3);
date_str=sprintf('%4.4i/%2.2i/%2.2i ',yr,mo,dm);

stt=sprintf('NEMO SSH %s',date_str);
%  if cntr==1
%    figure('Position',[1064 229 1421 1113]); clf;
%  end
clr1=[0 0.6 1];
clr2=[0 1 0.8];

figure(1); clf;
set(gcf,'Position',[1000 481 1383 835]);

POS = [0.05 0.08 0.44 0.85; ...
							0.53 0.08  0.44 0.85];
fnb=1;
c1=-0.5;
c2=0.5;
hhN=[];
xl1=-98;
xl2=-78.4;
yl1=14.;
yl2=31.;
clb=1;
ps1=POS(1,:);
sub_plot_ssh3(ps1,sshn,LONN,LATN,hhN,c1,c2,stt,INN,xl1,xl2,yl1,yl2,clb);

% ==============
% HYCOM-TSIS
% ==============
ps1=POS(2,:);
stt2=sprintf('NEMO+GLORYS intrp HYCOM SSH %s',date_str);
clb=0;

sub_plot_ssh3(ps1,sshh,LON,LAT,hhN,c1,c2,stt2,INH,xl1,xl2,yl1,yl2,clb);
%    plot(LCN(1).xx,LCN(1).yy,'-','Color',clr1,'Linewidth',1.5);
%    plot(LCH(1).xx,LCH(1).yy,'-','Color',clr2,'Linewidth',1.5);
text(-96.5,15,flnm,'Fontsize',14,'Color',[1 0.9 0],'Interpreter','none');

btx = 'check_archv_ssh.m'; 
bottom_text(btx,'pwd',1,'Position',[0.08 0.05 0.4 0.05]);





