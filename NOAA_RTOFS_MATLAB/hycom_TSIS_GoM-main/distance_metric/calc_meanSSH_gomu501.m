% Compute daily mean SSH from 3-hr SSH with tides
% 0.04 HYCOM GOMu reanalysis 50.1
% HYCOM + NCODA Gulf of Mexico 1/25° Reanalysis
% Tidal forcing is included!
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

clear all
close all

HR=0;
esim='GoMu501'; % 

% Years to calculate:
ys = 2009;
ye = ys;

f_mat = 1;  % = 0 - do not save,
            % = 1 - save matlab file with transport, etc.;
            % = 2 - load saved mat file and continue from where the last day is 
            %      in the mat file till the last day = dy2


rg=9806;  % convert pressure to depth, m
huge=1e10;

mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(ys,4)==0; mday(2)=29; end;
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
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
  end
end


% -------------------------
% Directories:
% -------------------------
ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/nexsan/archive/GOMu0.04_501/topo/';
btx = 'calc_volFlx_GoMreanalysis.m';


% -------------------------
% Get grid and bath:
% My bathymetry, Model bathymetry:
% -------------------------
% Read topography:
ftopo = sprintf('%sdepth_GOMu0.04_03i.nc',pthtopo);
HH  = -1*squeeze(nc_varget(ftopo,'depth'));
HH(isnan(HH))=100;
alat = nc_varget(ftopo,'Latitude');
elon = nc_varget(ftopo,'Longitude');
[mm,nn]=size(HH);
HH(isnan(HH))=100;


SSHM = struct;
mold = -1;
irc  = 0;
nrc  = size(YRPLT,1);
for ip=1:nrc
  iyr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%    icyc=YRPLT(ip,3);
  dJ1=datenum(iyr,1,1);
  dnmb=dJ1+iday-1;
  DV=datevec(dnmb);

  yr   = DV(1);
  mo   = DV(2);
  mday = DV(3);

%keyboard
  fprintf('Reading: %s; esim=%s\n',datestr(dnmb),esim);

% Daily average fields to get rid of tides:
  ihr = 0;
  ssh = zeros(mm,nn);  
  tic;
  for HR=0:3:21
    pthi=sprintf('/nexsan/archive/GOMu0.04_501/data/netcdf/%4.4i/',iyr);
    fina=sprintf('%shycom_gomu_501_%i%2.2i%2.2i00_t%3.3i.nc',...
                pthi,DV(1:3),HR);

    if ~exist(fina,'file')
      fprintf('Missing %s\n',fina);
      continue
    end

    ihr = ihr+1;
    dmm = squeeze(nc_varget(fina,'surf_el'));
    ssh = ssh+dmm;
  end

  ssh=ssh./ihr;

  fprintf(' Processed 1 day: %6.4f min, min/max ssh = %8.4f/%8.4f m\n',...
            toc/60,nanmin(nanmin(ssh)),nanmax(nanmax(ssh)));
  fprintf('--- \n');

  if mold<=0; mold=mo; end;
  if mold~=mo 
    fmatout = sprintf('%sssh_daily_GOMu501_%4.4i%2.2i.mat',pthmat,yr,mold);
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'SSHM'); 
    irc = 0;
    SSHM = struct;
    mold = mo;
  end
  irc = irc+1;
  SSHM(irc).YR     = yr;
  SSHM(irc).Mo     = mo;
  SSHM(irc).day    = mday;
  SSHM(irc).nrec   = ihr;
  SSHM(irc).ssh_mn = ssh; 

end

fmatout = sprintf('%sssh_daily_GOMu501_%4.4i%2.2i.mat',pthmat,yr,mold);
fprintf('Saving %s\n',fmatout);
save(fmatout,'SSHM');

sprintf(' ============   All Done    =============\n');
  




