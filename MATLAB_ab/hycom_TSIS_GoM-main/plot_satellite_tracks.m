% Plot ssh tracks
% for specific dates
% in the GoM
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

pthdat = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/obs/aviso/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';

dnmb1 = datenum(2009,12,01);
dnmb2 = datenum(2009,12,09);
dj1   = datenum(2009,1,1);

TM = [dnmb1:dnmb2];
nt = length(TM);

SAT{1}='en';
SAT{2}='j1n';
SAT{3}='j2';


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
TRCK(1).ST(1).Name = 'Envisat';
TRCK(2).ST(2).Name = 'Jason1';
TRCK(3).ST(3).Name = 'Jason2';
for id=1:nt
  dnmb=TM(id);
  dv=datevec(dnmb);

  for ist=1:3
    sat=SAT{ist};

    fnm=sprintf('%sdt_global_%s_phy_vfec_l3_%i%2.2i%2.2i.nc',pthdat,sat,dv(1:3));
    fprintf('Reading %s\n',fnm);

    lat = nc_varget(fnm,'latitude');
    lon = nc_varget(fnm,'longitude');
    J=find(lon>180);
    lon(J)=lon(J)-360;
    time = nc_varget(fnm,'time');

    I=find(lat>lt1 & lat<lt2 & lon>ln1 & lon<ln2);
 
    TRCK(id).ST(ist).X = lon(I);
    TRCK(id).ST(ist).Y = lat(I);
    TRCK(id).TM        = dnmb;

  end
end

% Plot tracks:
CLR = [0 0.5 0.9; ...
       0.9 0.4 0; ...
       0   0.9 0.4];


for id=1:9
  dnmb = TRCK(id).TM;
  dv=datevec(dnmb);
  fprintf('Plotting %i/%i/%i\n',dv(1:3));

  figure(id); clf;
  axes('Position',[0.08 0.2 0.85 0.72]);
  hold on;
  contour(LON,LAT,HH,[0 0],'color','k','linewidth',2);

  for ist=1:3
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
  text(-97,30.5,stl,'Fontsize',14);

  axes('Position',[0.09 0.26 0.2 0.15]);
  hold on;
  plot([0.02 0.2],[0.1 0.1],'-','Linewidth',2,'Color',CLR(1,:));
  text(0.3,0.1,'Envisat');
  plot([0.02 0.2],[0.2 0.2],'-','Linewidth',2,'Color',CLR(2,:));
  text(0.3,0.2,'Jason1');
  plot([0.02 0.2],[0.3 0.3],'-','Linewidth',2,'Color',CLR(3,:));
  text(0.3,0.3,'Jason2');
  set(gca,'xtick',[],...
          'ytick',[],...
          'xlim',[0 0.6],...
          'ylim',[0 0.4]);


  btx='plot_satellite_tracks.m';
  bottom_text(btx,'pwd',1);


end
    





