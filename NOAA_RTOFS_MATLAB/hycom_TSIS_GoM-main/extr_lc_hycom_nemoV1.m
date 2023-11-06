% Hindcast experiments:
% Compare HYCOM data assimilative run 2011- June 2012 vs
% NEMO simulations
% Observations are created from NEMO fields
%
% OSSE's experiments: initial fields in May 2011 and Jan 2012
% from HYCOM data assimilative runs with NEMO fields 
% different experiments - different set of synthetic obs. from NEMO being assimilated
%
%
% Analysis of LC, steps:
% 1) extract LC contour extr_lc_hycom_nemoV1.m
% 2) same for nemo LC contour:  extr_lc_ssh_nemo.m if needed
% 3) calculate MHD: distance_metric/mhd_osse_hindcasts_hycom_nemoV1.m
% 4) Plot results: distance_metrics/

% Extract and save LC contours from HYCOM-TSIS and
% nemo simulations with new hindcasts and free run 
% specify individually which run need to extract
%
%  NEMO is extracted in extr_lc_hycom_nemo-V0.m 
%
% Compare LC contours from the HYCOM_TSIS hindcast
% assimilating NEMO 1/100 free running fields
% and NEMO LC contours
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
EXON(2:9) = 1; % select expt to be extracted,  #2 - ssh ???
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

btx = 'extr_lc_hycom_nemoV1.m';

ii=0;
EXPT = struct;
% Experiments:
% freerun
ii=ii+1;
EXPT(ii).Name = 'FreeRun';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/freerun/';

% No SSS or SST analysis:
% full field SSH + no pies  
ii=ii+1;
EXPT(ii).Name = '2DSSH noPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_full_nopies_newtsis/';
% AVISO tracks SSH + no pies 
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathsSSH noPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_nopies_newtsis/';
% “one track” SSH + no pies
ii=ii+1;
EXPT(ii).Name = '1swathSSH noPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_onetrack_nopies_newtsis/';
% AVISO tracks SSH + ugos pies (no SSS and no SST)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH ugosPIES noSSS noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_ugos_nosss_nosst_newtsis/';

% No SSS analysis but with SST analysis:
%No SSH track + pies distributed all over the GOM domain (1/30 points)
ii=ii+1;
EXPT(ii).Name = 'noSSH allGoMPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_nosla_gompies_nosss_newtsis/';
%AVISO tracks SSH + ugos pies (small area distribution)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH ugosPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_ugos_nosss_newtsis/';
% AVISO tracks SSH + extd pies (bigger area distribution)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH extendedPIES noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_extd_nosss_newtsis/';
% No SSH track + pies distributed all over the GOM domain  (1/60 points)
ii=ii+1;
EXPT(ii).Name = 'NoSSH allGoMPIES60 noSSS SST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_nosla_gompies60_nosss_newtsis/';

fhnd = 'hycom_tsis_expts.mat';
%%%save(fhnd,'EXPT');  now prepared in hindcast_info.m

load(fhnd);

Nruns = ii;

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



cntr=0;
% Read in HYCOM ssh from requested experiments:
Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;  
  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fprintf('%s %s\n',nmexp,fmatout);

  clear LCXY
  LCXY.Name = nmexp;
  LCXY.Pthdata = pthd1;

  for ii=1:nrc
    tic;

    yr   = YPLT(ii,1);
    mo   = YPLT(ii,2);
    dm   = YPLT(ii,3);
    dyr  = YPLT(ii,4);
    dnmb = YPLT(ii,5);
    iday = dyr;

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
  % Derive LC contour:
  % 
    dmm=ssh;
    dmm(INH==0)=nan; 
    bsgn = -1;
    f_stop = 0;
    LCH1 = identify_LC(LON,LAT,dmm,Bisol,'bsgn',bsgn,'err_stop',f_stop);
%keyboard

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

    fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);

    if mod(ii,30)==0 & f_mat>0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'LCXY','LCE');
    end

  end;

  fprintf('Finished, Saving %s\n',fmatout);
  save(fmatout,'LCXY','LCE');

end


