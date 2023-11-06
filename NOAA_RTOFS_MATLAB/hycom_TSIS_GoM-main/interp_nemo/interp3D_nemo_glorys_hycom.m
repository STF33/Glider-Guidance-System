% Wrapper script
% 2D field only - that do not require vertical
% interpolation onto NEMO vertical layers
%
% Interpolation of 3D fields - T, S, U
% Interpolation consists of several steps:
% 1) Interpolate NEMO field onto HYCOM-TSIS
% 2) Interpolate GLORYS onto HYCOM-TSIS
% 3) Combine NEMO + GLORYS on HYCOM-TSIS
% 
% Note all interpolation weights, points are precomputed 
% NEMO - HYCOM takes long to locate points
% 
% If relax HYCOM fields need to be created: 
% write interpolated fields into climatology HYCOM *a. *b files:
% see: write_hycom_clim.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

f_write = 1;
fldnm = 'temp';
%fldnm = 'saln';
%fldnm = 'uvel';
%fldnm = 'vvel';

Nlrs = 75;     % depth layers in NEMO

if f_write < 0
  fprintf(' \n\n !!!!  Data are not saved !!!!  f_write=%i\n\n',f_write);
end

dnmb = datenum(2011,7,9);  % interpolation date
DV = datevec(dnmb);

fprintf('Interpolating: %s %s\n\n',fldnm,datestr(dnmb));

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';

% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;


[DX,DY]=sub_dx_dy(LON,LAT);
%[XM,YM]=meshgrid([1:nn],[1:mm]);

% Find NEMO bounding indices:
fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
fprintf('Loading %s\n',fnemo_indx);
NMI = load(fnemo_indx);
Inm = NMI.IndxHYCOM;   % indices covered by NEMO
[JH,IH] = ind2sub(size(HH),Inm);
jbot = min(JH);
jtop = max(JH);
ilft = min(IH);
irht = max(IH);

% 
% Load NEMO grid:
% Get NEMO grid:
f_get_grid=0;
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
if f_get_grid==1
  [LONN,LATN,ZZN] = sub_get_NEMO_grid(dnmb);
  LONN = double(LONN);
  LATN = double(LATN);
  ZZN = double(ZZN);
  fprintf('Saving grid %s\n',fgrd);
  save(fgrd,'LONN','LATN','ZZN');
else
  fprintf('Loading NEMO grid %s\n',fgrd);
  load(fgrd);
end


fprintf('%s interpolating into HYCOM-TSIS %s\n',fldnm,datestr(dnmb));

AmskN = [];  % T/S mask to fill bottom cells
AmskG = []; 
for iz0 = 1:Nlrs
%keyboard
  tic;
  zz0 = ZZN(iz0);  % NEMO depth
  fprintf('  %s iz=%i Z=%6.2f\n',fldnm,iz0,zz0);

%  if ZZN(iz0)<-250; keyboard; end;

% (1) interpolate NEMO to HYCOM
  [Anemo, AmskN] = sub_nemo2tsis_3d(dnmb,pthnemo,pthtopo,pthdata,...
                           LAT,LON,HH,DX,DY,fldnm,iz0,NMI,...
                           LONN,LATN,ZZN,AmskN);

% (2) interpolate GLORYS to HYCOM on NEMO depth layer
  [Aglrs, AmskG] = sub_glorys2tsis_3d(dnmb,pthnemo,pthdata,pthglorys,LAT,LON,HH,...
                             fldnm,zz0,AmskG);
%  if iz0>72; keyboard; end;
% 
% Merge two fields:


%Lmsk       = HH*0;
%Lmsk(HH<0) = 1;
%Lmsk(Inm)  = 0;

  Ahycom = Aglrs;
  Ahycom(Inm) = Anemo(Inm);
  A0 = Ahycom;

%
% Along NEMO bndry, there is discontinuity/jump 
% smooth it by applying Gaussian
% smooth only ourside of the NEMO domain
% fitting a spline along norm direction from the bndry
  nij = 7;  % total # of points around the bndry to interpolate
  din = 3;  % how many points in the domain
  Ahycom = sub_smooth_bndry(nij,din,Ahycom,jbot,jtop,ilft,irht);
% Filter
% di points, Rflt - Gaussian "radius"
  Rflt = 1;
  dout = nij;
  Ahycom = sub_filter_bndry(Rflt,Ahycom,jbot,jtop,ilft,irht,din,dout);

%
% Check discontinuity at the NEMO boundary:
  btx = 'interp3D_nemo_glorys_hycom.m';
  f_chck=0;
  if f_chck==1
    dltE=HH*nan;
    for ik=irht-20:irht+20
      d1=Ahycom(:,ik-1);
      d2=Ahycom(:,ik+1);
      dlt=d2-d1;
      dltE(jbot:jtop,ik)=dlt(jbot:jtop);
    end;

    for jk=jbot-20:jbot+20
      d1=Ahycom(jk-1,:);
      d2=Ahycom(jk+1,:);
      dlt=d2-d1;
      dltE(jk,ilft:irht)=dlt(ilft:irht);
    end;
      
    for jk=jtop-20:jtop+20
      d1=Ahycom(jk-1,:);
      d2=Ahycom(jk+1,:);
      dlt=d2-d1;
      dltE(jk,ilft:irht)=dlt(ilft:irht);
    end;

    c1=10;
    c2=30;

    if strncmp(fldnm,'saln',4);
      c1=30;
      c2=37;
    elseif strncmp(fldnm,'uvel',4) | strncmp(fldnm,'vvel',4)
      c1=-0.4;
      c2=0.4;
    end

    figure(2); clf;
    pcolor(Ahycom); shading flat;
    hold on;
    caxis([c1 c2]);
    axis('equal');
    
    plot([ilft irht],[jbot jbot],'k');
    plot([ilft ilft],[jbot jtop],'k');
    plot([irht irht],[jbot jtop],'k');
    plot([ilft irht],[jtop jtop],'k');
    colorbar('SouthOutside');
    title(sprintf('%s NEMO+GLORYS on HYCOM-TSIS grid',fldnm));
    bottom_text(btx,'pwd',1);

    di=30;
    figure(3); clf;
    hold on;
    contour(HH,[0 0],'-','Color',[0.6 0.6 0.6]);
    pcolor(dltE); shading flat;
    caxis([-0.02 0.02]);
    colorbar

    axis('equal');
    set(gca,'xlim',[400 irht+2*di]);
    set(gca,'ylim',[jbot-2*di jtop+2*di]);
    title('Jump across the NEMO bndry');

  %
  % Pick East North sections
    ik = 450; % East bndry
    jk = 605; % North bndry
    figure(5); clf;
    axes('Position',[0.09 0.53 0.84 0.38])
    hold on;
    dmm=Ahycom(ik,:);
    yl1=min(dmm(irht-di:irht+di));
    yl2=max(dmm(irht-di:irht+di));

    plot(dmm,'.-');

    if exist('A0','var'),
      dmm0=A0(ik,:);
      plot(dmm0,'-');
    end

    plot([irht irht],[min(dmm) max(dmm)],'r--');
    set(gca,'xlim',[irht-di irht+di]);
    set(gca,'ylim',[yl1 yl2]);
    title('East');

  % North
    axes('Position',[0.09 0.1 0.84 0.38])
    hold on;
    dmm=Ahycom(:,jk);
    yl1=min(dmm(jtop-di:jtop+di));
    yl2=max(dmm(jtop-di:jtop+di));

    plot(dmm,'.-');

    if exist('A0','var'),
      dmm0=A0(:,jk);
      plot(dmm0,'-');
    end

    plot([jtop jtop],[min(dmm) max(dmm)],'r--');
    set(gca,'xlim',[jtop-di jtop+di]);
    set(gca,'ylim',[yl1 yl2]);
    title('North');

    bottom_text(btx,'pwd',1);

    keyboard

  end


  % Fixed record length 
  [mm,nn]= size(HH);
  JDM = mm;
  IDM = nn;
  IJDM=IDM*JDM;
  npad=4096-mod(IJDM,4096);
  toto=ones(npad,1);

  if f_write==1
  % Fill land
    F = sub_fill_land(Ahycom);  
    F = F';
    F = reshape(F,IJDM,1);

    if isempty(fid)  
      fout = sprintf('%s%s_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                   pthoutp,fldnm,DV(3),DV(2),DV(1));
      fid  = fopen(fout,'w','ieee-be');
    end

    fwrite(fid,F,'float32','ieee-be');
    fwrite(fid,toto,'float32','ieee-be');
    fprintf('Written files: %s\n',fout);
    
  end

  fprintf(' --------  1 record %8.4f min\n\n',toc/60);
end;   % depth levels

if f_write==1
  fclose(fid);
end

%
% Update links to *a, and create *b climatology files that are needed
% for relax files (interpolated onto HYCOM hybrid layers T&S)
% This replaces the code write_hycom_clim.m
if f_write==1 & (strncmp(fldnm,'temp',4) | strncmp(fldnm,'saln',4))
  sgm = 'sig2';  % pressure reference

  flnma = sprintf('%s_zlev_hycom_%s_%2.2i%2.2i%4.4i.a',fldnm,sgm,DV(3),DV(2),DV(1));
  fouta = sprintf('%s%s',pthoutp,flnma); 
  ib = max(strfind(fout,'/'));
  fnm_new = fout(ib+1:end);
  fprintf('Creating link %s ---> %s\n',fnm_new,flnma);
  D1 = pwd;

  s1=sprintf('cd %s',pthoutp(1:end-1));
  s2=sprintf('rm -f %s',flnma);
  s3=sprintf('ln -sf %s %s',fnm_new,flnma);
  s4=sprintf('cd %s',D1);

%  s1='who';
%  [stat,cout] = system(s1);
  eval(s1);
  system(s2);
  system(s3);
  eval(s4);

% Create *b file
  [JDM,IDM] = size(HH);
  IJDM = IDM*JDM;
  sub_create_climb(pthoutp,fldnm,flnma,dnmb,IDM,JDM,ZZN);

end;











