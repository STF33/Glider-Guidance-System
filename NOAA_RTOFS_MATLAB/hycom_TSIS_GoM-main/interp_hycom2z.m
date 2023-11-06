% For LC contour finding: interpolate T or S fields to 
% Z depth Z0
%
% Extract LC and LCEs for calculating MHD
% Using T fields (demeaned) at depth Z0, contour T anomaly
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
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/interp2grid
startup;

close all
clear

f_mat = 1; % save mat; =2 - load saved and finish missing dates
Z0 = -200;  % depth


EXON = zeros(9,1);
EXON(6) = 1; % select expt to be extracted,  #2 - ssh ???
%ixx = 3;  % hindcast/free run # expt



pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';

% HYCOM-TSIS hindcast experiments:
fhnd = 'hycom_tsis_expts.mat';
load(fhnd);
%Nruns = length(EXPT);
Nruns = 9; % <10 OSSE expts, for Predictability expts >=10 use *Prdct* codes

for ii=1:Nruns
  if EXON(ii)==0
    fprintf('%i : OFF    %s \n',ii,EXPT(ii).Name);
  else
    fprintf('%i : ON ---> %s \n',ii,EXPT(ii).Name);
  end
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
Iocn = find(INH==1 & HH<Z0);
clear XM YM


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

XX=LON(jts1:jts2,its1:its2);
YY=LAT(jts1:jts2,its1:its2);
HHs=HH(jts1:jts2,its1:its2);
Hnas=HHs;
[msnsb]=size(HHs);

dz=50;
ZMnas = [Z0+dz:-dz:Z0-dz];
iz0 = find(ZMnas==Z0);

knas=length(ZMnas);
%[mnas,nnas]=size(nlon);
mnas = m;
nnas = n;

% Create Land/bottom mask indices:
for kk=1:knas
		I=find(Hnas>=ZMnas(kk));
		ILAND(kk).I=I;
end

% 
TZH.Info = 'Interpolated T to fixed z0';
TZH.Z0   = Z0;
TZH.HH   = HHs;
TZH.LON  = XX;
TZH.LAT  = YY;
TM = [];

% Read in HYCOM ssh from requested experiments:
Iexpt = find(EXON==1);
for jj=1:length(Iexpt);
  ixx = Iexpt(jj);
  nmexp = EXPT(ixx).Name;
  pthd1 = EXPT(ixx).path;
%  fmatout = sprintf('%shycom_LCcontour_%2.2i.mat',pthmat,ixx);
  fmatout = sprintf('%shycom_t2Z%4.4i_hindcast%2.2i.mat',pthmat,abs(Z0),ixx);

  fprintf('%s %s\n',nmexp,fmatout);

  TZH.Name = nmexp;
  TZH.Pthdata = pthd1;

  if f_mat==2
    fprintf('\n\n !!! Continue from last saved record, loading %s\n',fmatout);
    load(fmatout);
    TM = TZH.TM;
  end


  irec = 0;
  for ii=1:nrc

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

    if f_mat==2
      if dnmb<=TM(end); 
        idt = find(TM==dnmb);       
        if ~isempty(idt) 
          fprintf('irec=%i, Date %s found, skipping\n',idt,datestr(dnmb));
        else
          fprintf('======  Missing Date %s =====\n',datestr(dnmb));
        end
        continue;
      end
      irec = length(TM);
      fprintf('Continue interpolating, irec=%i\n',irec);
    end
      
    fprintf('Expt = %i %s\n',ixx,datestr(dnmb));


    fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthd1,yr,iday);
    finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthd1,yr,iday);
    fin=fina;

    ie = exist(fin,'file');

    if ~ie
      fprintf('Missing: %s\n',fin);
      continue;
    end
%keyboard
    tic;
    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F/rg;
    F(F>1e10)=0;
    dH=squeeze(F); % note this is not U,V thickness, center grid
        % there are thck-u & thck-v 
    ih=find(dH<1e-1);
    dHs=dH(:,jts1:jts2,its1:its2);

    [ZM1,ZZ1] = sub_zz_zm(fina,finb,HH,'thknss',dH);
    ZM1(isnan(ZM1))=-100;
    ZZ=ZZ1(:,jts1:jts2,its1:its2);
    ZM=ZM1(:,jts1:jts2,its1:its2);

%keyboard
    fprintf('Reading %s\n',fina);
    fld = 'temp';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    T=F;
%    T=F(:,jts1:jts2,its1:its2);

    if ~exist('IFL','var');
      IFL=[];
      fifl = sprintf('%sIFL.mat',pthmat);
      load(fifl);
    end

% fill nans in the surface layer  
    t1=squeeze(T(1,:,:));
%    [t1,IFL]=sub_fill_land(t1,IFL);
    t1=sub_fill_land(t1);
%   fifl = sprintf('%sIFL.mat',pthmat);
%    save(fifl,'IFL');
    T(1,:,:)=t1;
    T0=T;
% Subsample domain:
    T=T0(:,jts1:jts2,its1:its2);

% Fill nans in the subsurface layers (bottom):
    T = sub_fill_3z(T);

% Interpolate onto NAS grid:    
% Interpolate into z levels
% as a 2D vert sections
    [msb,nsb]=size(XX);
				Tzi=[];
				for ii=1:nsb
						if mod(ii,100)==0
								fprintf('%s interp to z, done %5.2f%% ...\n',fld,ii/nsb*100);
						end

						y1=YY(:,ii);
						ZMh=squeeze(ZM(:,:,ii));
						F=squeeze(T(:,:,ii));
						Hs = HHs(:,ii);
						Hs(Hs>=0)=-1;

						Fi=sub_interp2z_2D(Hs,F,ZMnas,ZMh);
      Tzi(:,:,ii)=Fi;

    end

% Put land back
				for kk=1:knas
						il=ILAND(kk).I;
						Tzi(kk,il)=nan;
				end
    Tz0 = squeeze(Tzi(iz0,:,:));

% Check:
%    T=T0(:,jts1:jts2,its1:its2);
%    aa = squeeze(T(13,:,:));

    irec=irec+1;
    TZH.TM(irec) = dnmb;
    TZH.Tz(irec,:,:) = Tz0;

    fprintf('min/max Tintrp= %8.2f %8.2f\n',min(min(Tz0)),max(max(Tz0)));    
    fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);

    if mod(irec,15)==0 & f_mat>0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'TZH');
    end

  end

  fprintf('---- End ----   Saving %s\n',fmatout);
  save(fmatout,'TZH');

end

