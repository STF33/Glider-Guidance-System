% Extract LC and LCEs for calculating MHD
%
% NEMO data
%
% NEMO data:
%ncdump -h https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc
%fin='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/2011/GOLFO108-NAS00_1d_20110101_20110131_grid_T.nc';
%ncid = netcdf.open(fin);
%vid  = netcdf.inqVarID(ncid,'ssh');
%ssh  = netcdf.getVar(ncid,vid,'single');
%[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vid)
%varname = netcdf.inqVar(ncid,vid)
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 2; % save mat; =2 - load saved and finish missing dates
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

YPLT=[];
cc=0;
for iy=2011:2012
%for iy=2010:2010
  for dd=1:366
%for iy=2011:2011
%  for dd=223:223
    if iy==2011 & dd==1; continue; end;
    if iy==2011 & dd==366; continue; end;
%    if iy==2012 & dd>182,
%      break;
%    end
    dnmb=datenum(iy,1,1)+dd-1;
    dv=datevec(dnmb);
    cc=cc+1;
    YPLT(cc,1)=iy;
    YPLT(cc,2)=dv(2);
    YPLT(cc,3)=dv(3);
    YPLT(cc,4)=dd;
    YPLT(cc,5)=dnmb;
  end
end

nrc=cc;


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


fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';

fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);
if iy==2010,
  fmatout = sprintf('%sNEMO_LCcontour%i.mat',pthmat,iy);
end

cntr=0;
LONN=[];
LATN=[];
INN=[];

irc_last=0;
dlast = 0;
if f_mat==2
  fprintf('Loading saved %s\n',fmatout);
  load(fmatout);
  irc_last = length(LCXY.TM);
  dlast    = LCXY.TM(end);
  cntr     = irc_last;
  fprintf('Last saved record # %i, %s\n\n',irc_last,datestr(dlast));

end

for ii=1:nrc
  tic;

  yr   = YPLT(ii,1);
  mo   = YPLT(ii,2);
  dm   = YPLT(ii,3);
  dyr  = YPLT(ii,4);
  dnmb = YPLT(ii,5);
  iday = dyr;

  if dnmb<=dlast,
    fprintf(' Record exists, skipping %s\n',datestr(dnmb));
    continue;
  end

  if f_mat == 2; f_mat=1; end;

  DV = datevec(dnmb);
  dnmb1=datenum(yr,mo,1);
  dnmb2=dnmb1+32;
  v2=datevec(dnmb2);
  dnmb2=datenum(v2(1),v2(2),1);
  d2=dnmb2-datenum(yr,mo,1);

  fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
  fprintf('Reading NEMO: %s\n',fnemo);

  fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

  if isempty(LONN)
    fmesh=sprintf('%smesh_mask.nc',fpwd);
    dmm = ncread(fmesh,'nav_lon');
    LONN = dmm';
    dmm = squeeze(ncread(fmesh,'nav_lat'));
    LATN = dmm';

    [mm,nn] = size(LONN);

    [XM,YM]=meshgrid([1:nn],[1:mm]);
    INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
    clear XM YM

  end


  enm = squeeze(ncread(fin,'ssh',[1 1 dm],[Inf Inf 1]));
  enm = enm';
  I=find(enm==0);
  enm(I)=nan;

  % Subtract spatial mean ssh
  dmm=enm;
  dmm(INN==0)=nan;
  sshM=nanmean(nanmean(dmm));
  enm = enm-sshM;

  Bisol=0.17;
  dmm=enm;
  dmm(INN==0)=nan;
  LCN = identify_LC(LONN,LATN,dmm,Bisol);

% Nemo
  cntr=cntr+1;
  LCXY.TM(cntr)    = dnmb;
  LCXY.XY(cntr).X  = LCN(1).xx;
  LCXY.XY(cntr).Y  = LCN(1).yy;
% Save LCEs as well: 
  LCE(1).TM(cntr)  = dnmb;
  lcc=length(LCN);
  LCE(1).NumbLCE(cntr) = lcc;
  LCE(1).XY(cntr).X=[];
  LCE(1).XY(cntr).Y=[];
  if lcc>1
    for ilc=2:lcc
      LCE(ilc-1).XY(cntr).X = LCN(ilc).xx;
      LCE(ilc-1).XY(cntr).Y = LCN(ilc).yy;
    end
  end

  fprintf('Processed 1 rec, %6.4f min\n',toc/60);

  if mod(ii,30)==0 & f_mat==1
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'LCXY','LCE');
  end

end
if f_mat==1
  fprintf('Saving NEMO %s\n',fmatout);
  save(fmatout,'LCXY','LCE');
end



return
