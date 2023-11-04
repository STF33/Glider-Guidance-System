% Wrapper script
% 2D field only - that do not require vertical
% interpolation onto NEMO vertical layers
%
% Interpolation of 2D fields - SSH
% Interpolation consists of several steps:
% 1) Interpolate NEMO field onto HYCOM-TSIS
% 2) Interpolate GLORYS onto HYCOM-TSIS
% 3) Combine NEMO + GLORYS on HYCOM-TSIS
% 
% Note all interpolation weights, points are precomputed 
% NEMO - HYCOM takes long to locate points
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_write = 1;
fldnm = 'ssh';

dnmb = datenum(2011,7,9);  % interpolation date - initial date of the f/cast
DV = datevec(dnmb);

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

fprintf('SSH interpolated into HYCOM-TSIS %s\n',datestr(dnmb));

% (1) interpolate NEMO to HYCOM
sshNi = sub_nemo2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,LAT,LON,HH,DX,DY);

% (2) interpolate GLORYS to HYCOM
sshGi = sub_glorys2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,pthglorys,LAT,LON,HH);

% 
% Merge two fields:
fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
fprintf('Loading %s\n',fnemo_indx);

% Find NEMO bounding indices:
NMI = load(fnemo_indx);
Inm = NMI.IndxHYCOM;   % indices covered by NEMO
[JH,IH] = ind2sub(size(HH),Inm);
jbot = min(JH);
jtop = max(JH);
ilft = min(IH);
irht = max(IH);

%Lmsk       = HH*0;
%Lmsk(HH<0) = 1;
%Lmsk(Inm)  = 0;

sshH = sshGi;
sshH(Inm) = sshNi(Inm);
sshH0=sshH;

%
% Along NEMO bndry, there is discontinuity/jump 
% smooth it by applying Gaussian
% smooth only ourside of the NEMO domain
% fitting a spline along norm direction from the bndry
nij = 7;  % total # of points around the bndry to interpolate
din = 3;  % how many points in the domain
sshH = sub_smooth_bndry(nij,din,sshH,jbot,jtop,ilft,irht);
% Filter
% di points, Rflt - Gaussian "radius"
Rflt = 1;
dout = nij;
sshH = sub_filter_bndry(Rflt,sshH,jbot,jtop,ilft,irht,din,dout);

%
% Check discontinuity at the NEMO boundary:
btx = 'interp2D_nemo_glorys_hycom.m';
f_chck=0;
if f_chck==1
  dltE=HH*nan;
  for ik=irht-20:irht+20
    d1=sshH(:,ik-1);
    d2=sshH(:,ik+1);
    dlt=d2-d1;
    dltE(jbot:jtop,ik)=dlt(jbot:jtop);
  end;

  for jk=jbot-20:jbot+20
    d1=sshH(jk-1,:);
    d2=sshH(jk+1,:);
    dlt=d2-d1;
    dltE(jk,ilft:irht)=dlt(ilft:irht);
  end;
    
  for jk=jtop-20:jtop+20
    d1=sshH(jk-1,:);
    d2=sshH(jk+1,:);
    dlt=d2-d1;
    dltE(jk,ilft:irht)=dlt(ilft:irht);
  end;

  figure(2); clf;
  pcolor(sshH); shading flat;
  hold on;
  caxis([-0.5 0.5]);
  axis('equal');
  
  plot([ilft irht],[jbot jbot],'k');
  plot([ilft ilft],[jbot jtop],'k');
  plot([irht irht],[jbot jtop],'k');
  plot([ilft irht],[jtop jtop],'k');
  colorbar('SouthOutside');
  title('SSH NEMO+GLORYS on HYCOM-TSIS grid');
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
  dmm=sshH(ik,:);
  yl1=min(dmm(irht-di:irht+di));
  yl2=max(dmm(irht-di:irht+di));

  plot(dmm,'.-');

  if exist('sshH0','var'),
    dmm0=sshH0(ik,:);
    plot(dmm0,'-');
  end

  plot([irht irht],[-0.3 0.5],'r--');
  set(gca,'xlim',[irht-di irht+di]);
  set(gca,'ylim',[yl1 yl2]);
  title('East');

% North
  axes('Position',[0.09 0.1 0.84 0.38])
  hold on;
  dmm=sshH(:,jk);
  yl1=min(dmm(jtop-di:jtop+di));
  yl2=max(dmm(jtop-di:jtop+di));

  plot(dmm,'.-');

  if exist('sshH0','var'),
    dmm0=sshH0(:,jk);
    plot(dmm0,'-');
  end

  plot([jtop jtop],[-0.3 0.5],'r--');
  set(gca,'xlim',[jtop-di jtop+di]);
  set(gca,'ylim',[yl1 yl2]);
  title('North');

  bottom_text(btx,'pwd',1);

end


% Fixed record length 
[mm,nn]= size(HH);
JDM = mm;
IDM = nn;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

if f_write==1
  clear sshGi sshNi sshH0
% Fill land
  F = sub_fill_land(sshH);  
  F = F';
  F = reshape(F,IJDM,1);

  fout = sprintf('%s%s_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                 pthoutp,fldnm,DV(3),DV(2),DV(1));
  fid  = fopen(fout,'w','ieee-be');
  fwrite(fid,F,'float32','ieee-be');
  fwrite(fid,toto,'float32','ieee-be');
  fclose(fid);
  fprintf('Written files: %s\n',fout);
  
end



