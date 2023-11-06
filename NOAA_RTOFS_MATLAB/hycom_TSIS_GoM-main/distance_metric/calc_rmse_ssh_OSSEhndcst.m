% Calculate RMSE for ssh
% sqrt (sum[(a-<a>)^2]/n) 
%  OSSE hindcast and NEMO 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
% Set flags for extracting experiments:
EXON = zeros(9,1);
EXON(5) = 1; % select expt to be extracted,  #2 - ssh ???
f_nemo = 0;    % obsolete - use extr_lc_ssh_nemo.m
Bisol = 0.17;  % ssh contour

if f_nemo>0
  fprintf('NEMO is now extracted in extr_lc_ssh_nemo.m !!! flag ignored\n\n');
  f_nemo=0;
end

% Hindcasts:
pthd1  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_new/';
pthd2  = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_extd_new/';  % extended PIES arrays
pthd12 = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/hindcast_ugos_obs/';  % all satellites included
pthd3  = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2011_GLfreerun/'; % free run
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

%fmatout = sprintf('%sLC_coord_osse_hycomV1.mat',pthmat);

btx = 'calc_rmse_ssh.m';


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


YPLT=[];
cc=0;
for iy=2011:2012
  for dd=1:365
    if iy==2011 & dd==1; continue; end;
    if iy==2012 & dd>182,
      break;
    end
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


%
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


% NEMO
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
LONN=[];
LATN=[];
INN=[];

% Find closest nemo index
f_indx = 0;
findxmat = sprintf('%snemo_hycom_indx.mat',pthmat);
if f_indx==1
  dnmb = datenum(2011,03,01);
  DV = datevec(dnmb);
  yr=DV(1);
  mo=DV(2);
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

  II=find(INH==1);
  np=length(II);
  H2N.I=INH*nan;
  H2N.J=INH*nan;
  H2N.linIndx=II*nan;
  mnD=0;
  fprintf('Searching NEMO indices\n');
  for ik=1:np
    if mod(ik,round(np/20))==0
      fprintf('mean dx=%6.2fkm  %5.1f%% done ...\n',mnD/ik*1e-3,ik/np*100);
    end
    i0=II(ik); 
    x0=LON(i0);
    y0=LAT(i0);
    D=distance_spheric_coord(LATN,LONN,y0,x0);
    mind=min(min(D));
    if mind>4000, continue; end; % outside NEMO
    [jn,in]=find(D==mind);
    H2N.I(i0)=in;
    H2N.J(i0)=jn;
    H2N.linIndx(i0)=sub2ind(size(D),jn,in);

    mnD=mnD+mind;
  end

  fprintf('Saving %s\n',findxmat);
  save(findxmat,'H2N');

else
  fprintf('Loading %s\n',findxmat);
  load(findxmat);
end;


%keyboard

% Read in HYCOM ssh from requested experiments:
Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;  
  fmatout = sprintf('%shycom_nemo_RMSE_hcst%2.2i.mat',pthmat,ixx);
  fprintf('%s %s\n',nmexp,fmatout);

  clear RMSE
  RMSE.Name = nmexp;
  RMSE.Pthdata = pthd1;


  Ihn=H2N.I;
  Jhn=H2N.J;
  II0=find(~isnan(Ihn));
  Ihycom=II0;
  Inemo=H2N.linIndx(II0);

  sumR = INH*0;
  Iout = find(INH==0);
  sumR(Iout)=nan;
  cntr=0;

  if f_mat==2
    fprintf('Loading %s\n',fmatout);
    load(fmatout);
    TMs=RMSE.TM;
%keyboard
  end

  for ii=1:nrc
    tic;

    yr   = YPLT(ii,1);
    mo   = YPLT(ii,2);
    dm   = YPLT(ii,3);
    dyr  = YPLT(ii,4);
    dnmb = YPLT(ii,5);
    iday = dyr;

    if f_mat==2
      if dnmb<=TMs(end);
        cntr=cntr+1;
        fprintf('Already processed, skipping %s\n',datestr(dnmb));
        if dnmb==TMs(end)
          sumR=RMSE.rmse.^2*cntr;
        end
        continue
      end
    end

    dnmb1=datenum(yr,mo,1);
    dnmb2=dnmb1+32;
    v2=datevec(dnmb2);
    dnmb2=datenum(v2(1),v2(2),1);
    d2=dnmb2-datenum(yr,mo,1);

    sday=sprintf('%3.3i',iday);
    hr=0;

    
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

%  
%       NEMO
% 
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


  % RMSE:
    cntr=cntr+1;
    dmm=(ssh(Ihycom)-enm(Inemo)).^2;
    sumR(Ihycom)=sumR(Ihycom)+dmm;

    mxR=max(max(sqrt(sumR/cntr)));
    fprintf('Processed 1 rec, %6.4f min, hcst# %i, maxRMSE=%8.4fcm\n\n',toc/60,ixx,mxR);

    RMSE.TM(cntr)=dnmb;
    if mod(ii,30)==0 & f_mat>0
      RMSE.rmse = sqrt(sumR/cntr);
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'RMSE');
    end

  end;
  RMSE.rmse = sqrt(sumR/cntr);
  
  fprintf('Finished, Saving %s\n',fmatout);
  save(fmatout,'RMSE');

end


