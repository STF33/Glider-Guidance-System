% Plot ssh tracks
% for specific dates
% in the GoM - IAS domain
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

s_mat = 1; 

%pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/obs/aviso/';
pthdat  = '/nexsan/people/abozec/TSIS/data/aviso/';  % <-- all data here, gzip
%pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/data/aviso/';  % copied from Alex's dir
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/obs/aviso/';

dnmb1 = datenum(2011,1,1);
dnmb2 = datenum(2012,7,1);

TM = [dnmb1:dnmb2];
nt = length(TM);

% For 2009;
%SAT{1}='en';
%SAT{2}='j1n';
%SAT{3}='j2';

% For 2011-2012
% c2 - Cryosat
%  enn - Envisat
% j1g - Jason 1 Geodetic Phase Global Ocean
% j1n - Jason 1 interleaved 
% j2 - Jason 2
%
SAT{1}='c2';
SAT{2}='enn';
SAT{3}='j1g';
SAT{4}='j1n';
SAT{5}='j2';
nst = length(SAT);


% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

ln1=min(min(LON));
ln2=max(max(LON));
lt1=min(min(LAT));
lt2=max(max(LAT));

TRCK = struct;
TRCK(1).ST(1).Name = 'Cryosat';
TRCK(1).ST(2).Name = 'Envisat';
TRCK(1).ST(3).Name = 'Jsn1Geod';
TRCK(1).ST(4).Name = 'Jsn1Intrlv';
TRCK(1).ST(5).Name = 'Jason2';


for id=1:nt
  dnmb=TM(id);
  dv=datevec(dnmb);

  TRCK(id).TimeMat = dnmb;
  TRCK(id).TimeStr = datestr(dnmb);

  for ist=1:nst
    sat=SAT{ist};
    fnm=sprintf('%s%s/dt_global_%s_phy_vfec_l3_%i%2.2i%2.2i.nc',pthdat,sat,sat,dv(1:3));

    if ~exist(fnm,'file') 
      fprintf('%s does not exist ...\n',fnm);
      continue;
    end

    fprintf('Reading %s\n',fnm);

    lat = nc_varget(fnm,'latitude');
    lon = nc_varget(fnm,'longitude');
    J=find(lon>180);
    lon(J)=lon(J)-360;
    time = nc_varget(fnm,'time');

    I=find(lat>lt1 & lat<lt2 & lon>ln1 & lon<ln2);
 
    TRCK(id).ST(ist).X = lon(I);
    TRCK(id).ST(ist).Y = lat(I);

  end
end

%keyboard

fmatout = sprintf('%saviso_L3_sattracks_2011-2012.mat',pthmat);
if s_mat==1
  fprintf('\n\n =====>  Saving %s\n',fmatout);
  save(fmatout,'TRCK');
end

f_chck=0;

if f_chck==1
  load(fmatout);
  dnmb = datenum(2011,09,01);
  

% Plot tracks:
  CLR = [0 0.5 0.9; ...
       0.9 0.4 0; ...
       0   1 0.7; ...
       0.8 0 1; ...
       0.4 0.6 0.2];

  nl = length(TRCK);

  for id=1:nl
    if TRCK(id).TimeMat == dnmb; break; end;
  end
  dv=datevec(dnmb);
  fprintf('Plotting %i/%i/%i\n',dv(1:3));

  figure(1); clf;
  axes('Position',[0.08 0.2 0.85 0.72]);
  hold on;
  contour(LON,LAT,HH,[0 0],'color','k','linewidth',2);

  for ist=1:nst
    x=TRCK(id).ST(ist).X;
    y=TRCK(id).ST(ist).Y;
  
    clr = CLR(ist,:);
    plot(x,y,'.','Markersize',12,'Color',clr);
  end

  axis('equal');
  set(gca,'tickdir','out',...
          'xlim',[ln1 ln2],...
          'ylim',[lt1 31.9]);
  stl=sprintf('%i/%2.2i/%2.2i',dv(1:3));
%  title(stl,'Fontsize',12);
  text(-97,30.55,stl,'Fontsize',14);

  axes('Position',[0.1 0.26 0.25 0.18]);
  hold on;
  for ist=1:nst
    snm = TRCK(1).ST(ist).Name;
    plot([0.02 0.15],[ist/10 ist/10],'-','Linewidth',2,'Color',CLR(ist,:));
    text(0.3,ist/10,snm);
  end
  set(gca,'xtick',[],...
          'ytick',[],...
          'xlim',[0 0.55],...
          'ylim',[0 0.55]);


  btx='extract_satellite_tracks.m';
  bottom_text(btx,'pwd',1);


end
    





