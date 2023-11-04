% Extract LC to calculate MHD
% OSE hindcasts 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

ys = 2009;
ye = 2009;
esim = 'noPIES';
%esim = 'PIES';
esim = 'PIESv2';  % updated TSIS with fixed assim problems of T/S prof
Z0 = -200;

f_mat = 1; % =1 save and overide existing mat; =2 - load saved and finish missing dates

if f_mat==0
  fprintf('\nOutput is not saved !!! \n\n');
  fprintf('Check f_mat=%i\n',f_mat);
end

%
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthout  = '/Net/kronos/ddmitry/hycom/TSIS/datafcst1b/';
pthfrcst = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_osse_fcst/';

% HYCOM:
rg=9806;  % convert pressure to depth, m
huge=1e20;
Bisol = 0.17;  % ssh contour

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
Iocn = find(INH==1 & HH<Z0);
Nocn=length(Iocn);
clear XM YM

Imean = INH;
Imean(HH>Z0)=0;


% 
% Subsample to a smaller domain:
xnas1 = min(LON(Iocn));
xnas2 = max(LON(Iocn));
ynas1 = min(LAT(Iocn));
ynas2 = max(LAT(Iocn));
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);


YRPLT=[];
cc=0;
for iyr=ys:ye
  id1=1;
  id2=365;
  ic=mod(iyr,4);
  if ic==0 & id2==365, id2=366; end;
  if iyr>ys
    id1=1;
  end
  for iday=id1:id2
    cc=cc+1;
    jd1=datenum(iyr,1,1);
    dnmb=jd1+iday-1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
    YRPLT(cc,3)=dnmb;
    DV(cc,:) = datevec(dnmb);
  end
end
nrc=size(YRPLT,1);

fmatout = sprintf('%shycom_LCcontour_OSEhindcast_%s_%2.2i-%2.2i.mat',pthmat,esim,ys,ye);

%  pthi2=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_noPIES/',iyr);

fprintf('\n\n %s\n',fmatout);
fprintf(' OSE hindcast: %i- %i \n',ys,ye);
fprintf(' Output saved: %s\n',fmatout);

dmean=1;
cntr=0;
for ii=1:nrc
	tic;
	yr  = DV(ii,1);
	mo  = DV(ii,2);
	dm  = DV(ii,3);
	dnmb= YRPLT(ii,3);
	iday= YRPLT(ii,2);
	HR  = 0;
  pthi=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%4.4i_%s/',yr,esim);
  if strncmp(esim,'PIESv2',6)
    pthi='/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_newtsis/gofs30_withpies/';
  end

	fina = sprintf('%sarchv.%i_%3.3i_%2.2i.a',pthi,yr,iday,HR);
	finb = sprintf('%sarchv.%i_%3.3i_%2.2i.b',pthi,yr,iday,HR);

	if ~exist(fina,'file')
		fprintf('Missing %s\n',fina);
		continue
	end

	ssh = sub_getSSH_hycom(fina,finb,Imean,dmean);

%
% Derive LC contour:
% 
	dmm=ssh;
	dmm(INH==0)=nan;
	bsgn = -1;
	f_stop = 0;
	LCH1 = identify_LC(LON,LAT,dmm,Bisol,'err_stop',f_stop);

	cntr=cntr+1;

% HYCOM-TSIS
	LCXY.TM(cntr)    = dnmb;
	LCXY.XY(cntr).X  = LCH1(1).xx;
	LCXY.XY(cntr).Y  = LCH1(1).yy;
% Save LCEs as well:
	LCE(1).TM(cntr)  = dnmb;
	lcc=length(LCH1);
	LCE(1).NumbLCE(cntr) = lcc;
	LCE(1).XY(cntr).X=[];
	LCE(1).XY(cntr).Y=[];
	if lcc>1
		for ilc=2:lcc
			LCE(ilc-1).XY(cntr).X = LCH1(ilc).xx;
			LCE(ilc-1).XY(cntr).Y = LCH1(ilc).yy;
		end
	end

	fprintf('cntr=%i, Processed 1 rec, %6.4f min\n\n',cntr,toc/60);

	btx = 'extr_lc_OSEhindcast.m';
	f_chck = 0;
	if ii==20 & f_chck==1
		figure(10); clf;
		pcolor(LON,LAT,dmm); shading flat;
		hold on;
		contour(LON,LAT,dmm,[Bisol Bisol],'k');

		axis('equal');
		set(gca,'xlim',[-98 -80],...
						'ylim',[17 31]);

		xlc = LCXY.XY(cntr).X;
		ylc = LCXY.XY(cntr).Y;
		dstr = datestr(dnmb);
		plot(xlc,ylc,'r.');

		nlce = length(LCE);
		for jj=1:nlce
		 xlce = LCE(jj).XY(cntr).X;
		 ylce = LCE(jj).XY(cntr).Y;
		 if isempty(xlce), continue; end;

		 plot(xlce,ylce,'m.');
		end
		stl = sprintf('%s %s',esim,dstr);
		title(stl);
		bottom_text(btx,'pwd',1);

		keyboard
	end

	if mod(ii,20)==0 & f_mat>0
			fprintf('Saving %s\n',fmatout);
			save(fmatout,'LCXY','LCE');
	end

end

fprintf('======   Finished,  Saving %s\n',fmatout);
save(fmatout,'LCXY','LCE');


