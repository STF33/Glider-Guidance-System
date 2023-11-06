% check SSH and LC/LCE contours from control and f/cast
% used for MHD calculation
% SSH of the f/cast is plotted
%
function sub_plot_ssh_LCLCE_chck(iFcst,itime,irun,dnmb,XYfcst,XYcntr,FCST,stl)

%IEXX = zeros(9,1);  % hindcast free run # expt
%IEXX(2:9) = 1;
%ixx = 10;  % hindcast/free run # expt
%iFcst = ixx;
% Specify forecast run:
%nhnd  = iFcst; %initial conditions, 10 - interpolation
%itime = 1;  % =1 2011, =2 2012
%irun  = 1;  % f/cast # within the time period - 7-day shifted
%
% Plot for:
%irc1 = 53;
%irc2 = 53;

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthmat2 = '/Net/kronos/ddmitry/hycom/TSIS/datafcst2/';
pthfrm = '/Net/kronos/ddmitry/hycom/TSIS/FIG/frames_fcst2/';


% HYCOM-TSIS hindcast experiments:
%fhnd = 'hycom_tsis_expts.mat';
%load(fhnd);
%Nruns = length(EXPT);

%
% Get LC/LCE contours from the f/cast
%fmatout = sprintf('%shycom_LCcontour_fcst%2.2i-%2.2i%2.2i.mat',pthmat2,nhnd,itime,irun);
%fprintf('Loading %s\n',fmatout);
%load(fmatout);
%TM  = LCXY.TM;
%nrec = length(TM);
%LCHY=LCXY;
%LCEH=LCE;

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

% 
% HYCOM output f/cast:
%if iFcst<10
%  FCST = sub_fcst_info(iFcst);
%else
%  FCST = sub_fcstPrdct_info(iFcst);
%end
%
pthhycom = FCST.TIME0(itime).RUN(irun).pthbin; % output dir

% GoM region HYCOM:
% similar to NEMO GOM domain
% for getting same endpoints of the contours
GOM=[366   489
   476   531
   542   560
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



cmp=flipud(colormap_cold(360));

figure(10); 
set(gcf,'Position',[1364         502         944         840]);


DV   = datevec(dnmb);
yr   = DV(1);
mo   = DV(2);
dm   = DV(3);
iday = dnmb-datenum(yr,1,1)+1;


% HYCOM   
% Get HYCOM ssh:
fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthhycom,yr,iday);
finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthhycom,yr,iday);
fprintf('Reading %s\n',fina);

huge  = 2e20;
rg    = 9806;

fld = 'srfhgt';
[F,nn,mm,ll] = read_hycom(fina,finb,fld);
F(F>huge)=nan;
ssh=squeeze(F)./(1e-3*rg);  % ssh m
dmm=ssh;
dmm(INH==0)=nan;
%  dmm(HH>-200)=nan;
sshM=nanmean(nanmean(dmm));
ssh=ssh-sshM;

%
% Get contours from NEMO and HYCOM hindcast:
%nmexp = EXPT(ixx).Name;
%pthd1 = EXPT(ixx).path;

%inm = find(TMN==dnmb);
%dm1 = TMN(inm);

% Control run contours:
xcntr = XYcntr(:,1); % 
ycntr = XYcntr(:,2);;


% F/cast contours:
xfcst = XYfcst(:,1);
yfcst = XYfcst(:,2);

%
% PLOTTING
clf;
axes('Position',[0.08 0.2 0.8 0.7]);
pcolor(LON,LAT,ssh); shading flat;
%stl=sprintf('ssh, MHD=%6.4d HYCOM fcst %i  %2.2i/%2.2i/%4.4i',mhd,iFcst,DV(3),DV(2),DV(1));
colormap(cmp);
caxis([-0.5 0.5]);
hold;

plot(xcntr,ycntr,'.','Color',[0 0 0]);
plot(xfcst,yfcst,'.','Color',[1 0 0]);

contour(LON,LAT,HH,[0 0],'k');

plot([-97.5 -96.5],[30.2 30.2],'r-','Linewidth',2.0);
text(-96.4,30.2,sprintf('F/cast '));
plot([-97.5 -96.5],[29.8 29.8],'k-','Linewidth',2.0);
text(-96.4,29.6,sprintf('Control'));


axis('equal');
set(gca,'xlim',[-97.8 -80.7],...
	'ylim',[18.5 31],...
	'xtick',[-98:2:-82],...
	'ytick',[18:2:32],...
	'tickdir','out',...
	'Fontsize',12);
title(stl);

hb=colorbar('SouthOutside');
set(hb,'Position',[0.25 0.15 0.5 0.015],...
 'Ticks',[-1:0.25:1],...
 'Fontsize',13);





